//
// Created by etoga on 4/11/23.
//

#ifndef DYNAMICCONVEXHULL_AVLTREE_H
#define DYNAMICCONVEXHULL_AVLTREE_H

#include <algorithm>
#include <iostream>
#include <vector>
#include <vector>

using uint = unsigned int;

template<class T, typename Comparator = std::less<T>>
class AVLTree{
    Comparator less;

public:
    struct Node{
        Node *left = nullptr, *right = nullptr, *par = nullptr;
        uint size = 1;
        int height = 1;
        T val, min, max;

        Node() = default;

        explicit Node(T v) : val(v),min(v),max(v) {};

        // Deletion of entire subtree
        ~Node(){
            delete left;
            delete right;
        }

        inline
        void makeLeftChild(Node* p){
            if(!p) return;
            p->par = this;
            left = p;
        }

        inline
        void makeRightChild(Node* p){
            if(!p) return;
            p->par = this;
            right = p;
        }
    };

    Node* root = nullptr;

    // Default constructor
    AVLTree() = default;
    
    // Destructor
    ~AVLTree(){
        delete root;
    }
    
    // Deep copy constructor - creates a complete independent copy of the tree
    AVLTree(const AVLTree& other) : root(deepCopyNode(other.root)) {}
    
    // Deep copy assignment operator
    AVLTree& operator=(const AVLTree& other) {
        if (this != &other) {
            delete root;
            root = deepCopyNode(other.root);
        }
        return *this;
    }
    
private:
    // Helper to recursively deep copy a node and its subtree
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
    
    // Move constructor
    AVLTree(AVLTree&& other) noexcept : root(other.root) {
        other.root = nullptr;
    }
    
    // Move assignment operator
    AVLTree& operator=(AVLTree&& other) noexcept {
        if (this != &other) {
            delete root;
            root = other.root;
            other.root = nullptr;
        }
        return *this;
    }

    void insert(T val){
        Node* parent = nullptr;
        Node* current = root;
        while(current){
            onVisit(current);
            parent = current;
            if(isLeaf(current)) break;
            else if(less(val,current->right->min)) current = current->left;
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
        if(!current || current->val != val) return; // Val not found
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

    void join(T k, AVLTree* r){
        Node* tl = root;
        Node* tr = r->root;
        Node* x = new Node(k);
        if(height(tl) > height(tr)+1) joinRight(root,x, r->root);
        else if(height(tr) > height(tl)+1) joinLeft(root,x,r->root);
        else {
           x->makeLeftChild(tl);
           x->makeRightChild(tr);
           updateData(x);
        }
        while(x->par) x = x->par;
        root = x;
        r->root = nullptr;
    }

    void join(Node* k, AVLTree* r) {
        Node *tl = root;
        Node *tr = r->root;
        Node *x = k;
        
        // In leaf-oriented trees, internal nodes MUST have exactly two children.
        // If one subtree is empty, we cannot create a valid internal node.
        // Instead, delete the pivot and return the non-empty subtree directly.
        
        if (!tl && !tr) {
            // Both empty - the pivot is unused, delete it
            // Result is an empty tree
            delete x;
            root = nullptr;
            r->root = nullptr;
            return;
        }
        if (!tl) {
            // Only right tree exists - delete pivot, use right tree
            delete x;
            root = tr;
            if (root) root->par = nullptr;
            r->root = nullptr;
            return;
        }
        if (!tr) {
            // Only left tree exists - delete pivot, keep left tree
            delete x;
            // root is already tl, just ensure no parent pointer
            if (root) root->par = nullptr;
            r->root = nullptr;
            return;
        }
        
        // Both trees are non-empty - create the join with pivot
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
        join(temp,r);

    }

    // Split s.t. l < k and r > k
    // After split: 'this' contains nodes < k, 'r' contains nodes > k
    // If k exists (as a leaf), it is deleted
    void split(T k, AVLTree* r){
        if(!root) {
            if(r) r->root = nullptr;
            return;
        }
        
        // Use recursive helper
        Node* rightRoot = nullptr;
        root = splitRecursive(root, k, rightRoot);
        
        // Clean up parent pointers for new roots
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
    // Recursively update all nodes in a subtree, bottom-up
    void updateAllNodes(Node* node) {
        if (!node) return;
        // First update children (bottom-up)
        updateAllNodes(node->left);
        updateAllNodes(node->right);
        // Then update this node
        updateData(node);
    }
    // Recursive split helper that returns the left tree root and sets rightRoot to the right tree root
    // The node 'node' is the current subtree to split
    // Returns: root of the < k portion
    // Sets rightRoot to: root of the > k portion
    Node* splitRecursive(Node* node, const T& k, Node*& rightRoot) {
        if(!node) {
            rightRoot = nullptr;
            return nullptr;
        }
        
        if(isLeaf(node)) {
            // Base case: leaf node
            if(node->val == k) {
                // Key matches - delete the leaf
                delete node;
                rightRoot = nullptr;
                return nullptr;
            } else if(less(node->val, k)) {
                // Leaf value < k, so it goes to left tree
                node->par = nullptr;
                rightRoot = nullptr;
                return node;
            } else {
                // Leaf value > k, so it goes to right tree
                node->par = nullptr;
                rightRoot = node;
                return nullptr;
            }
        }
        
        // Optimization checking bounds (min/max)
        // If node entirely < k, returns {node, nullptr}
        if (less(node->max, k)) {
             rightRoot = nullptr;
             node->par = nullptr;
             // Internal structure preserved, bridges valid internally.
             return node;
        }
        // If node entirely >= k, returns {nullptr, node}
        if (!less(node->min, k)) {
             rightRoot = node;
             node->par = nullptr;
             // Internal structure preserved, bridges valid internally.
             return nullptr;
        }

        // Internal node that must be split
        // 1. Recursive split
        Node* leftLeft = nullptr;
        Node* leftRight = nullptr;
        // Optimization: if k is in right subtree, entire left child goes left.
        // But we just use recursion for simplicity first, or check bounds.
        // We know node->right->min is the separator? No, node->min and node->max are.
        
        // Standard recursive split
        Node* LL = splitRecursive(node->left, k, leftRight);  // leftRight is LR
        Node* RL = splitRecursive(node->right, k, rightRoot); // rightRoot is RR, RL is returned
        
        // 2. Capture pivot value before deleting node
        T pivot = node->val;
        
        // 3. Delete old internal node (structure dissolved)
        node->left = nullptr;
        node->right = nullptr;
        delete node;
        
        // 4. Reassemble using join
        // New Left Tree: join(LL, RL)
        Node* newLeft = joinNodes(LL, RL, pivot);
        
        // New Right Tree: join(LR, RR)
        rightRoot = joinNodes(leftRight, rightRoot, pivot); // LR, RR
        
        return newLeft;
    }

    // Helper to join two subtrees L and R using a dummy pivot
    // Guarantees proper recomputation of bridges via updateData
    Node* joinNodes(Node* L, Node* R, const T& dummyVal) {
        if (!L && !R) return nullptr;
        if (!L) return R;
        if (!R) return L;
        
        // Both exist: join them into a new internal node
        // Create new internal node using dummyVal as pivot
        Node* P = new Node(dummyVal);
        
        // Use standard join logic based on heights
        if (height(L) > height(R) + 1) {
            joinRight(L, P, R);
        } else if (height(R) > height(L) + 1) {
            joinLeft(L, P, R);
        } else {
            P->makeLeftChild(L);
            P->makeRightChild(R);
            updateData(P);
        }
        
        // Retrace up to find new root
        while (P->par) P = P->par;
        return P;
    }

protected:

    void joinRight(Node* tl, Node* x, Node* tr){
        while(tl->right && height(tl) > height(tr) + 1) tl = tl->right;
        if(tl && tl->par) tl->par->makeRightChild(x);
        x->makeLeftChild(tl);
        x->makeRightChild(tr);
        retrace(x,true, true);
    }
    void joinLeft(Node* tl, Node* x, Node* tr){
        while(tr->left && height(tr) > height(tl) + 1) tr = tr->left;
        if(tr && tr->par) tr->par->makeLeftChild(x);
        x->makeLeftChild(tl);
        x->makeRightChild(tr);
        retrace(x,true, true);
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
            //if ((insertion && dif == 0) || (!insertion && (dif == -1 || dif == 1))) balance = false;
        }
    }

    inline
    void updateData(Node* x){
        x->size = size(x->left) + size(x->right) + isLeaf(x);
        x->height = std::max(height(x->left), height(x->right)) + 1;
        if(x->left) x->min = x->left->min;
        else x->min = x->val;
        if(x->right) x->max = x->right->max;
        else x->max = x->val;
        onUpdate(x);
    }

    virtual inline
    void onUpdate(Node* x){};

    virtual inline
    void onVisit(Node* x){};

    inline
    bool isLeaf(const Node* x){
        return x && !(x->left || x->right);
    }

    inline
    bool isLeftChild(const Node* x){
        return (x->par && x->par->left == x);
    }

    inline
    int height(const Node* x){
        if(x) return x->height;
        return 0;
    }

    inline
    uint size(const Node* x){
        if(x) return x->size;
        return 0;
    }

public:
    // Clear all nodes and reset the tree
    void clear() {
        delete root;
        root = nullptr;
    }

    // Build a balanced AVL tree from a sorted vector in O(n) time.
    // The vector MUST be sorted according to the Comparator.
    // This replaces any existing tree content.
    void buildFromSorted(const std::vector<T>& sorted) {
        clear();
        if (sorted.empty()) return;
        root = buildFromSortedRecursive(sorted, 0, sorted.size());
    }

protected:
    // Recursive helper that builds a subtree from sorted[start..end)
    // Returns the root of the subtree.
    // For a leaf-oriented AVL tree: internal nodes store "bridges" computed from children,
    // and leaves store the actual data.
    Node* buildFromSortedRecursive(const std::vector<T>& sorted, size_t start, size_t end) {
        if (start >= end) return nullptr;
        
        if (end - start == 1) {
            // Base case: single element becomes a leaf
            Node* leaf = new Node(sorted[start]);
            return leaf;
        }
        
        // Internal node: pick median point for balanced split
        size_t mid = start + (end - start) / 2;
        
        // Create an internal node (its val will be computed by onUpdate)
        Node* node = new Node(sorted[mid]); // Temporary value, will be updated
        
        // Recursively build left and right subtrees
        node->left = buildFromSortedRecursive(sorted, start, mid);
        node->right = buildFromSortedRecursive(sorted, mid, end);
        
        // Set parent pointers
        if (node->left) node->left->par = node;
        if (node->right) node->right->par = node;
        
        // Update node data (height, size, min, max, and call onUpdate for bridges)
        updateData(node);
        
        return node;
    }
};

#endif //DYNAMICCONVEXHULL_AVLTREE_H
