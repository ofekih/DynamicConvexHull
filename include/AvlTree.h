// Copyright 2026 DynamicConvexHull Authors
// SPDX-License-Identifier: MIT

/// @file AvlTree.h
/// @brief Leaf-oriented AVL tree with split and join operations.
///
/// This AVL tree stores values only at leaves. Internal nodes store aggregate
/// information (min, max) and support efficient split/join operations.
///
/// Time Complexities:
///   - insert/remove: O(log n)
///   - split/join: O(log n)
///   - BuildFromSorted: O(n)

#pragma once

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <vector>

namespace dch {

/// @brief Leaf-oriented AVL tree with split/join support.
/// @tparam T Value type stored in leaves.
/// @tparam Comparator Comparison functor for ordering.
template <class T, typename Comparator = std::less<T>>
class AVLTree {
 public:
  /// @brief Tree node storing value and aggregate information.
  struct Node {
    Node* left = nullptr;
    Node* right = nullptr;
    Node* par = nullptr;
    std::size_t size = 1;
    int height = 1;
    T val;
    T min;
    T max;

    Node() = default;
    explicit Node(T v) : val(v), min(v), max(v) {}

    ~Node() {
      delete left;
      delete right;
    }

    void MakeLeftChild(Node* child) {
      if (!child) return;
      child->par = this;
      left = child;
    }

    void MakeRightChild(Node* child) {
      if (!child) return;
      child->par = this;
      right = child;
    }
  };

  AVLTree() = default;

  ~AVLTree() { delete root_; }

  /// @brief Deep copy constructor.
  AVLTree(const AVLTree& other) : root_(DeepCopyNode(other.root_)) {}

  /// @brief Deep copy assignment.
  AVLTree& operator=(const AVLTree& other) {
    if (this != &other) {
      delete root_;
      root_ = DeepCopyNode(other.root_);
    }
    return *this;
  }

  /// @brief Move constructor.
  AVLTree(AVLTree&& other) noexcept : root_(other.root_) {
    other.root_ = nullptr;
  }

  /// @brief Move assignment.
  AVLTree& operator=(AVLTree&& other) noexcept {
    if (this != &other) {
      delete root_;
      root_ = other.root_;
      other.root_ = nullptr;
    }
    return *this;
  }

  [[nodiscard]] Node* GetRoot() const { return root_; }

  /// @brief Insert a value. O(log n)
  void Insert(T val) {
    Node* parent = nullptr;
    Node* current = root_;
    while (current) {
      OnVisit(current);
      parent = current;
      if (IsLeaf(current)) break;
      if (less_(val, current->right->min)) {
        current = current->left;
      } else {
        current = current->right;
      }
    }
    Node* inserted = new Node(val);
    inserted->par = parent;
    if (!parent) {
      root_ = inserted;
      return;
    }
    Node* sibling = new Node(parent->val);
    sibling->par = parent;
    if (less_(val, parent->val)) {
      parent->left = inserted;
      parent->right = sibling;
    } else {
      parent->right = inserted;
      parent->left = sibling;
    }

    Retrace(parent, /*insertion=*/true, /*early_exit=*/false);
  }

  /// @brief Remove a value. O(log n)
  void Remove(T val) {
    Node* parent;
    Node* current = root_;
    while (current) {
      OnVisit(current);
      parent = current;
      if (IsLeaf(current)) break;
      if (less_(val, current->right->min)) {
        current = current->left;
      } else {
        current = current->right;
      }
    }
    if (!current || current->val != val) return;
    if (current == root_) {
      delete current;
      root_ = nullptr;
      return;
    }
    current = parent;
    parent = parent->par;
    Node* sibling;
    if (IsLeftChild(current)) {
      sibling = parent->right;
    } else {
      sibling = parent->left;
    }
    sibling->par = parent->par;
    if (!parent->par) {
      root_ = sibling;
    } else if (IsLeftChild(parent)) {
      parent->par->left = sibling;
    } else {
      parent->par->right = sibling;
    }

    current->left = current->right = nullptr;
    delete current;
    parent->left = parent->right = nullptr;
    delete parent;

    OnVisit(sibling);
    Retrace(sibling, /*insertion=*/false, /*early_exit=*/false);
  }

  /// @brief Join with another tree using a pivot value.
  /// @param pivot Pivot value to use as separator.
  /// @param other Tree to join (becomes empty after operation).
  void Join(T pivot, AVLTree* other) {
    Node* tree_left = root_;
    Node* tree_right = other->root_;
    Node* x = new Node(pivot);
    if (Height(tree_left) > Height(tree_right) + 1) {
      JoinRight(root_, x, other->root_);
    } else if (Height(tree_right) > Height(tree_left) + 1) {
      JoinLeft(root_, x, other->root_);
    } else {
      x->MakeLeftChild(tree_left);
      x->MakeRightChild(tree_right);
      UpdateData(x);
    }
    while (x->par) x = x->par;
    root_ = x;
    other->root_ = nullptr;
  }

  /// @brief Join with another tree using an existing node as pivot.
  void Join(Node* pivot_node, AVLTree* other) {
    Node* tree_left = root_;
    Node* tree_right = other->root_;
    Node* x = pivot_node;

    if (!tree_left && !tree_right) {
      delete x;
      root_ = nullptr;
      other->root_ = nullptr;
      return;
    }
    if (!tree_left) {
      delete x;
      root_ = tree_right;
      if (root_) root_->par = nullptr;
      other->root_ = nullptr;
      return;
    }
    if (!tree_right) {
      delete x;
      if (root_) root_->par = nullptr;
      other->root_ = nullptr;
      return;
    }

    if (Height(tree_left) > Height(tree_right) + 1) {
      JoinRight(root_, x, other->root_);
    } else if (Height(tree_right) > Height(tree_left) + 1) {
      JoinLeft(root_, x, other->root_);
    } else {
      x->MakeLeftChild(tree_left);
      x->MakeRightChild(tree_right);
      UpdateData(x);
    }
    while (x->par) x = x->par;
    root_ = x;
    other->root_ = nullptr;
  }

  /// @brief Join with another tree (extracts pivot from this tree).
  void Join(AVLTree* other) {
    Node* current = root_;
    if (!current) {
      root_ = other->root_;
      other->root_ = nullptr;
      return;
    }
    if (!other->root_) return;
    while (current->right) current = current->right;
    T temp = current->val;
    Split(current->val, nullptr);
    Join(temp, other);
  }

  /// @brief Split tree at key.
  ///
  /// After split: 'this' contains nodes < key, 'other' contains nodes > key.
  /// If key exists as a leaf, it is deleted.
  ///
  /// @param key Split key.
  /// @param other Receiver for right portion (may be nullptr to discard).
  void Split(T key, AVLTree* other) {
    if (!root_) {
      if (other) other->root_ = nullptr;
      return;
    }

    Node* right_root = nullptr;
    root_ = SplitRecursive(root_, key, right_root);

    if (root_) {
      root_->par = nullptr;
    }

    if (other) {
      other->root_ = right_root;
      if (right_root) {
        right_root->par = nullptr;
      }
    }
  }

  /// @brief Clear all nodes.
  void Clear() {
    delete root_;
    root_ = nullptr;
  }

  /// @brief Build balanced tree from sorted values in O(n) time.
  /// @param sorted Values sorted by Comparator.
  void BuildFromSorted(const std::vector<T>& sorted) {
    Clear();
    if (sorted.empty()) return;
    root_ = BuildFromSortedRecursive(sorted, 0, sorted.size());
  }

 protected:
  Node* root_ = nullptr;
  Comparator less_;

  static Node* DeepCopyNode(const Node* node) {
    if (!node) return nullptr;
    Node* copy = new Node(node->val);
    copy->size = node->size;
    copy->height = node->height;
    copy->min = node->min;
    copy->max = node->max;
    copy->left = DeepCopyNode(node->left);
    copy->right = DeepCopyNode(node->right);
    if (copy->left) copy->left->par = copy;
    if (copy->right) copy->right->par = copy;
    return copy;
  }

  void UpdateAllNodes(Node* node) {
    if (!node) return;
    UpdateAllNodes(node->left);
    UpdateAllNodes(node->right);
    UpdateData(node);
  }

  Node* SplitRecursive(Node* node, const T& key, Node*& right_root) {
    if (!node) {
      right_root = nullptr;
      return nullptr;
    }

    if (IsLeaf(node)) {
      if (node->val == key) {
        delete node;
        right_root = nullptr;
        return nullptr;
      } else if (less_(node->val, key)) {
        node->par = nullptr;
        right_root = nullptr;
        return node;
      } else {
        node->par = nullptr;
        right_root = node;
        return nullptr;
      }
    }

    if (less_(node->max, key)) {
      right_root = nullptr;
      node->par = nullptr;
      return node;
    }
    if (!less_(node->min, key)) {
      right_root = node;
      node->par = nullptr;
      return nullptr;
    }

    Node* left_right = nullptr;
    Node* left_left = SplitRecursive(node->left, key, left_right);
    Node* right_left = SplitRecursive(node->right, key, right_root);

    T pivot = node->val;

    node->left = nullptr;
    node->right = nullptr;
    delete node;

    Node* new_left = JoinNodes(left_left, right_left, pivot);
    right_root = JoinNodes(left_right, right_root, pivot);

    return new_left;
  }

  Node* JoinNodes(Node* left, Node* right, const T& dummy_val) {
    if (!left && !right) return nullptr;
    if (!left) return right;
    if (!right) return left;

    Node* pivot = new Node(dummy_val);

    if (Height(left) > Height(right) + 1) {
      JoinRight(left, pivot, right);
    } else if (Height(right) > Height(left) + 1) {
      JoinLeft(left, pivot, right);
    } else {
      pivot->MakeLeftChild(left);
      pivot->MakeRightChild(right);
      UpdateData(pivot);
    }

    while (pivot->par) pivot = pivot->par;
    return pivot;
  }

  void JoinRight(Node* tree_left, Node* x, Node* tree_right) {
    while (tree_left->right && Height(tree_left) > Height(tree_right) + 1) {
      tree_left = tree_left->right;
    }
    if (tree_left && tree_left->par) tree_left->par->MakeRightChild(x);
    x->MakeLeftChild(tree_left);
    x->MakeRightChild(tree_right);
    Retrace(x, /*insertion=*/true, /*early_exit=*/true);
  }

  void JoinLeft(Node* tree_left, Node* x, Node* tree_right) {
    while (tree_right->left && Height(tree_right) > Height(tree_left) + 1) {
      tree_right = tree_right->left;
    }
    if (tree_right && tree_right->par) tree_right->par->MakeLeftChild(x);
    x->MakeLeftChild(tree_left);
    x->MakeRightChild(tree_right);
    Retrace(x, /*insertion=*/true, /*early_exit=*/true);
  }

  void RotateRight(Node* x) {
    Node* y = x->left;
    OnVisit(y);
    x->left = y->right;
    if (y->right) y->right->par = x;
    y->par = x->par;
    if (!x->par) {
      root_ = y;
    } else if (IsLeftChild(x)) {
      x->par->left = y;
    } else {
      x->par->right = y;
    }
    y->right = x;
    x->par = y;

    UpdateData(x);
    UpdateData(y);
  }

  void RotateLeft(Node* x) {
    Node* y = x->right;
    OnVisit(y);
    x->right = y->left;
    if (y->left) y->left->par = x;
    y->par = x->par;
    if (!x->par) {
      root_ = y;
    } else if (IsLeftChild(x)) {
      x->par->left = y;
    } else {
      x->par->right = y;
    }
    y->left = x;
    x->par = y;

    UpdateData(x);
    UpdateData(y);
  }

  void Rebalance(Node* x) {
    OnVisit(x);
    if (Height(x->right) > Height(x->left)) {
      if (Height(x->right->right) < Height(x->right->left)) {
        OnVisit(x->right);
        RotateRight(x->right);
      }
      RotateLeft(x);
    } else if (Height(x->right) < Height(x->left)) {
      if (Height(x->left->right) > Height(x->left->left)) {
        OnVisit(x->left);
        RotateLeft(x->left);
      }
      RotateRight(x);
    }
  }

  void Retrace(Node* x, const bool insertion, const bool early_exit) {
    bool balance = true;
    int diff;

    for (auto current = x; current; current = current->par) {
      UpdateData(current);
      if (!balance) {
        continue;
      }

      diff = Height(current->left) - Height(current->right);

      if (diff < -1 || diff > 1) {
        Rebalance(current);
        current = current->par;
        diff = Height(current->left) - Height(current->right);
        if (insertion || diff == -1 || diff == 1) {
          balance = false;
          continue;
        }
      }
    }
  }

  void UpdateData(Node* x) {
    x->size = Size(x->left) + Size(x->right) + (IsLeaf(x) ? 1 : 0);
    x->height = std::max(Height(x->left), Height(x->right)) + 1;
    if (x->left) {
      x->min = x->left->min;
    } else {
      x->min = x->val;
    }
    if (x->right) {
      x->max = x->right->max;
    } else {
      x->max = x->val;
    }
    OnUpdate(x);
  }

  virtual void OnUpdate(Node* x) {}
  virtual void OnVisit(Node* x) {}

  [[nodiscard]] bool IsLeaf(const Node* x) const {
    return x && !(x->left || x->right);
  }

  [[nodiscard]] static bool IsLeftChild(const Node* x) {
    return (x->par && x->par->left == x);
  }

  [[nodiscard]] static int Height(const Node* x) {
    if (x) return x->height;
    return 0;
  }

  [[nodiscard]] static std::size_t Size(const Node* x) {
    if (x) return x->size;
    return 0;
  }

  Node* BuildFromSortedRecursive(const std::vector<T>& sorted,
                                  std::size_t start, std::size_t end) {
    if (start >= end) return nullptr;

    if (end - start == 1) {
      return new Node(sorted[start]);
    }

    std::size_t mid = start + (end - start) / 2;

    Node* node = new Node(sorted[mid]);

    node->left = BuildFromSortedRecursive(sorted, start, mid);
    node->right = BuildFromSortedRecursive(sorted, mid, end);

    if (node->left) node->left->par = node;
    if (node->right) node->right->par = node;

    UpdateData(node);

    return node;
  }
};

}  // namespace dch
