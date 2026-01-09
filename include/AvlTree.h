/**
 * @file AvlTree.h
 * @brief Leaf-oriented AVL tree with split and join operations.
 * 
 * This AVL tree stores values only at leaves. Internal nodes store aggregate
 * information (min, max) and support efficient split/join operations.
 * 
 * Time Complexities:
 *   - insert/remove: O(log n)
 *   - split/join: O(log n)
 *   - buildFromSorted: O(n)
 */

#ifndef DYNAMICCONVEXHULL_AVLTREE_H
#define DYNAMICCONVEXHULL_AVLTREE_H

#include <algorithm>
#include <iostream>
#include <vector>

using uint = unsigned int;

/**
 * @brief Leaf-oriented AVL tree with split/join support.
 * @tparam T Value type stored in leaves.
 * @tparam Comparator Comparison functor for ordering.
 */
template<class T, typename Comparator = std::less<T>>
class AVLTree{
    Comparator less;

public:
    /**
     * @brief Tree node storing value and aggregate information.
     */
    struct Node{
        Node *left = nullptr, *right = nullptr, *par = nullptr;
        uint size = 1;
        int height = 1;
        T val, min, max;

        Node() = default;
        explicit Node(T v) : val(v), min(v), max(v) {}

        ~Node(){
            delete left;
            delete right;
        }

        inline void makeLeftChild(Node* p){
            if(!p) return;
            p->par = this;
            left = p;
        }

        inline void makeRightChild(Node* p){
            if(!p) return;
            p->par = this;
            right = p;
        }
    };

    Node* root = nullptr;

    AVLTree() = default;
    
    ~AVLTree(){
        delete root;
    }
    
    /** @brief Deep copy constructor. */
    AVLTree(const AVLTree& other) : root(deepCopyNode(other.root)) {}
    
    /** @brief Deep copy assignment. */
    AVLTree& operator=(const AVLTree& other) {
        if (this != &other) {
            delete root;
            root = deepCopyNode(other.root);
        }
        return *this;
    }
    
private:
    static Node* deepCopyNode(const Node* node) {
        if (!node) return nullptr;
        Node* copy = new Node(node->val);
        copy->size = node->size;
        copy->height = node->height;
        copy->min = node->min;
        copy->max = node->max;
        copy->left = deepCopyNode(node->left);
        copy->right = deepCopyNode(node->right);
        if (copy->left) copy->left->par = copy;
        if (copy->right) copy->right->par = copy;
        return copy;
    }
    
public:
    Node* getRoot() const { return root; }
    
    /** @brief Move constructor. */
    AVLTree(AVLTree&& other) noexcept : root(other.root) {
        other.root = nullptr;
    }
    
    /** @brief Move assignment. */
    AVLTree& operator=(AVLTree&& other) noexcept {
        if (this != &other) {
            delete root;
            root = other.root;
            other.root = nullptr;
        }
        return *this;
    }

    /** @brief Insert a value. O(log n) */
    void insert(T val){
        Node* parent = nullptr;
        Node* current = root;
        while(current){
            onVisit(current);
            parent = current;
            if(isLeaf(current)) break;
            else if(less(val, current->right->min)) current = current->left;
            else current = current->right;
        }
        Node* inserted = new Node(val);
        inserted->par = parent;
        if(!parent){
            root = inserted;
            return;
        }
        Node* sibling = new Node(parent->val);
        sibling->par = parent;
        if(less(val, parent->val)){
            parent->left = inserted;
            parent->right = sibling;
        } else {
            parent->right = inserted;
            parent->left = sibling;
        }

        retrace(parent, true, false);
    }

    /** @brief Remove a value. O(log n) */
    void remove(T val){
        Node* parent;
        Node* current = root;
        while(current){
            onVisit(current);
            parent = current;
            if(isLeaf(current)) break;
            else if(less(val, current->right->min)) current = current->left;
            else current = current->right;
        }
        if(!current || current->val != val) return;
        if(current == root){
            delete current;
            root = nullptr;
            return;
        }
        current = parent;
        parent = parent->par;
        Node* sibling;
        if(isLeftChild(current)) sibling = parent->right;
        else sibling = parent->left;
        sibling->par = parent->par;
        if(!parent->par) root = sibling;
        else if(isLeftChild(parent)) parent->par->left = sibling;
        else parent->par->right = sibling;

        current->left = current->right = nullptr;
        delete current;
        parent->left = parent->right = nullptr;
        delete parent;

        onVisit(sibling);
        retrace(sibling, false, false);
    }

    /**
     * @brief Join with another tree using a pivot value.
     * @param k Pivot value to use as separator.
     * @param r Tree to join (becomes empty after operation).
     */
    void join(T k, AVLTree* r){
        Node* tl = root;
        Node* tr = r->root;
        Node* x = new Node(k);
        if(height(tl) > height(tr)+1) joinRight(root, x, r->root);
        else if(height(tr) > height(tl)+1) joinLeft(root, x, r->root);
        else {
           x->makeLeftChild(tl);
           x->makeRightChild(tr);
           updateData(x);
        }
        while(x->par) x = x->par;
        root = x;
        r->root = nullptr;
    }

    /**
     * @brief Join with another tree using an existing node as pivot.
     */
    void join(Node* k, AVLTree* r) {
        Node *tl = root;
        Node *tr = r->root;
        Node *x = k;
        
        if (!tl && !tr) {
            delete x;
            root = nullptr;
            r->root = nullptr;
            return;
        }
        if (!tl) {
            delete x;
            root = tr;
            if (root) root->par = nullptr;
            r->root = nullptr;
            return;
        }
        if (!tr) {
            delete x;
            if (root) root->par = nullptr;
            r->root = nullptr;
            return;
        }
        
        if (height(tl) > height(tr) + 1) joinRight(root, x, r->root);
        else if (height(tr) > height(tl) + 1) joinLeft(root, x, r->root);
        else {
            x->makeLeftChild(tl);
            x->makeRightChild(tr);
            updateData(x);
        }
        while (x->par) x = x->par;
        root = x;
        r->root = nullptr;
    }

    /** @brief Join with another tree (extracts pivot from this tree). */
    void join(AVLTree* r){
        Node* current = root;
        if(!current){
            root = r->root;
            r->root = nullptr;
            return;
        } else if(!r->root) return;
        while(current->right) current = current->right;
        T temp = current->val;
        split(current->val, nullptr);
        join(temp, r);
    }

    /**
     * @brief Split tree at key k.
     * 
     * After split: 'this' contains nodes < k, 'r' contains nodes > k.
     * If k exists as a leaf, it is deleted.
     * 
     * @param k Split key.
     * @param r Receiver for right portion (may be nullptr to discard).
     */
    void split(T k, AVLTree* r){
        if(!root) {
            if(r) r->root = nullptr;
            return;
        }
        
        Node* rightRoot = nullptr;
        root = splitRecursive(root, k, rightRoot);
        
        if(root) {
            root->par = nullptr;
        }
        
        if(r) {
            r->root = rightRoot;
            if(rightRoot) {
                rightRoot->par = nullptr;
            }
        }
    }

protected:
    void updateAllNodes(Node* node) {
        if (!node) return;
        updateAllNodes(node->left);
        updateAllNodes(node->right);
        updateData(node);
    }

    Node* splitRecursive(Node* node, const T& k, Node*& rightRoot) {
        if(!node) {
            rightRoot = nullptr;
            return nullptr;
        }
        
        if(isLeaf(node)) {
            if(node->val == k) {
                delete node;
                rightRoot = nullptr;
                return nullptr;
            } else if(less(node->val, k)) {
                node->par = nullptr;
                rightRoot = nullptr;
                return node;
            } else {
                node->par = nullptr;
                rightRoot = node;
                return nullptr;
            }
        }
        
        if (less(node->max, k)) {
             rightRoot = nullptr;
             node->par = nullptr;
             return node;
        }
        if (!less(node->min, k)) {
             rightRoot = node;
             node->par = nullptr;
             return nullptr;
        }

        Node* leftRight = nullptr;
        Node* LL = splitRecursive(node->left, k, leftRight);
        Node* RL = splitRecursive(node->right, k, rightRoot);
        
        T pivot = node->val;
        
        node->left = nullptr;
        node->right = nullptr;
        delete node;
        
        Node* newLeft = joinNodes(LL, RL, pivot);
        rightRoot = joinNodes(leftRight, rightRoot, pivot);
        
        return newLeft;
    }

    Node* joinNodes(Node* L, Node* R, const T& dummyVal) {
        if (!L && !R) return nullptr;
        if (!L) return R;
        if (!R) return L;
        
        Node* P = new Node(dummyVal);
        
        if (height(L) > height(R) + 1) {
            joinRight(L, P, R);
        } else if (height(R) > height(L) + 1) {
            joinLeft(L, P, R);
        } else {
            P->makeLeftChild(L);
            P->makeRightChild(R);
            updateData(P);
        }
        
        while (P->par) P = P->par;
        return P;
    }

    void joinRight(Node* tl, Node* x, Node* tr){
        while(tl->right && height(tl) > height(tr) + 1) tl = tl->right;
        if(tl && tl->par) tl->par->makeRightChild(x);
        x->makeLeftChild(tl);
        x->makeRightChild(tr);
        retrace(x, true, true);
    }

    void joinLeft(Node* tl, Node* x, Node* tr){
        while(tr->left && height(tr) > height(tl) + 1) tr = tr->left;
        if(tr && tr->par) tr->par->makeLeftChild(x);
        x->makeLeftChild(tl);
        x->makeRightChild(tr);
        retrace(x, true, true);
    }

    void rotateR(Node* x){
        Node* y = x->left;
        onVisit(y);
        x->left = y->right;
        if(y->right) y->right->par = x;
        y->par = x->par;
        if(!x->par) root = y;
        else if(isLeftChild(x)) x->par->left = y;
        else x->par->right = y;
        y->right = x;
        x->par = y;

        updateData(x);
        updateData(y);
    }

    void rotateL(Node* x){
        Node* y = x->right;
        onVisit(y);
        x->right = y->left;
        if(y->left) y->left->par = x;
        y->par = x->par;
        if(!x->par) root = y;
        else if(isLeftChild(x)) x->par->left = y;
        else x->par->right = y;
        y->left = x;
        x->par = y;

        updateData(x);
        updateData(y);
    }

    void rebalance(Node* x){
        onVisit(x);
        if(height(x->right) > height(x->left)){
            if(height(x->right->right) < height(x->right->left)){
                onVisit(x->right);
                rotateR(x->right);
            }
            rotateL(x);
        } else if (height(x->right) < height(x->left)){
            if(height(x->left->right) > height(x->left->left)) {
                onVisit(x->left);
                rotateL(x->left);
            }
            rotateR(x);
        }
    }

    void retrace(Node* x, const bool insertion, const bool earlyExit){
        bool balance = true;
        int dif;

        for(auto current = x; current; current = current->par){
            updateData(current);
            if(!balance){
                continue;
            }

            dif = height(current->left) - height(current->right);

            if(dif < -1 || dif > 1){
                rebalance(current);
                current = current->par;
                dif = height(current->left) - height(current->right);
                if(insertion || dif == -1 || dif == 1){
                    balance = false;
                    continue;
                }
            }
        }
    }

    inline void updateData(Node* x){
        x->size = size(x->left) + size(x->right) + isLeaf(x);
        x->height = std::max(height(x->left), height(x->right)) + 1;
        if(x->left) x->min = x->left->min;
        else x->min = x->val;
        if(x->right) x->max = x->right->max;
        else x->max = x->val;
        onUpdate(x);
    }

    virtual inline void onUpdate(Node* x){}
    virtual inline void onVisit(Node* x){}

    inline bool isLeaf(const Node* x){
        return x && !(x->left || x->right);
    }

    inline bool isLeftChild(const Node* x){
        return (x->par && x->par->left == x);
    }

    inline int height(const Node* x){
        if(x) return x->height;
        return 0;
    }

    inline uint size(const Node* x){
        if(x) return x->size;
        return 0;
    }

public:
    /** @brief Clear all nodes. */
    void clear() {
        delete root;
        root = nullptr;
    }

    /**
     * @brief Build balanced tree from sorted values in O(n) time.
     * @param sorted Values sorted by Comparator.
     */
    void buildFromSorted(const std::vector<T>& sorted) {
        clear();
        if (sorted.empty()) return;
        root = buildFromSortedRecursive(sorted, 0, sorted.size());
    }

protected:
    Node* buildFromSortedRecursive(const std::vector<T>& sorted, size_t start, size_t end) {
        if (start >= end) return nullptr;
        
        if (end - start == 1) {
            Node* leaf = new Node(sorted[start]);
            return leaf;
        }
        
        size_t mid = start + (end - start) / 2;
        
        Node* node = new Node(sorted[mid]);
        
        node->left = buildFromSortedRecursive(sorted, start, mid);
        node->right = buildFromSortedRecursive(sorted, mid, end);
        
        if (node->left) node->left->par = node;
        if (node->right) node->right->par = node;
        
        updateData(node);
        
        return node;
    }
};

#endif //DYNAMICCONVEXHULL_AVLTREE_H
