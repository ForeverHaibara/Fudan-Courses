#include<iostream>
using namespace std;

#define RBNode Node<T>*
#define RED 1
#define BLACK 0

template<class T>
class Node{
public:
    RBNode parent;
    RBNode left;
    RBNode right;
    int s;
    T val;
    bool color;
    Node(RBNode _p=nullptr,RBNode _l=nullptr,RBNode _r=nullptr,int _s=0,T _v=0,bool _c=1):
    parent(_p),left(_l),right(_r),s(_s),val(_v),color(_c){}

    inline void update(){
        s = 1 + left->s + right->s;
    }
    inline bool isleft(){
        return this == parent->left;
    }
    inline bool isright(){
        return this == parent->right;
    }
    inline bool isnil(){
        return this == left;
    }
    ~Node(){
        if (!(left->isnil())){
            delete left;
            left = nullptr;
        }
        if (!(right->isnil())){
            delete right;
            right = nullptr;
        }
    }
};

template<class T>
class RBTree{
private:
    RBNode root;
    RBNode nil;
    void doprint(RBNode x);
    void doprintStruct(RBNode x);
    RBNode get_lower_bound(RBNode x, T v);
    RBNode get_upper_bound(RBNode x, T v);
    inline void left_rotate(RBNode x);
    inline void right_rotate(RBNode x);
    void FixDoubleRed(RBNode x);
    void FixDoubleBlack(RBNode x);

public:
    RBTree(){
        root = nil->parent = nil->left = nil->right = nil = new Node<T>();
        nil->s = nil->color = nil->val = 0;
    }

    inline void print();
    inline void printStruct();
    RBNode successor(RBNode x);
    RBNode predecessor(RBNode x);
    RBNode lower_bound(T x);
    RBNode upper_bound(T x);
    inline RBNode rfind(T x);
    inline RBNode findkth(int x);
    inline int get_rank(T x);
    void insert(T x);
    void remove(RBNode x);

    ~RBTree(){
        // C++ delete空指针不报错
        delete root;    
        root = nullptr;
        delete nil;
        nil = nullptr;
    }
};

template<class T>
void RBTree<T>::doprint(RBNode x){
    if (x == nil) return;
    doprint(x->left);
    cout << x->val << ' ';
    doprint(x->right);
}

template<class T>
inline void RBTree<T>::print(){
    cout << "Tree: ";
    doprint(root);
    cout << endl;
}

template<class T>
void RBTree<T>::doprintStruct(RBNode x){
    if (x == nil) return;
    cout << (x->color==RED?'R':'B') << x->val << '[';
    doprintStruct(x->left);
    cout << ','; 
    doprintStruct(x->right); 
    cout << ']';
}

template<class T>
inline void RBTree<T>:: printStruct(){
    cout << "Tree: ";
    doprintStruct(root);
    cout << endl;
}

template <class T>
RBNode RBTree<T>:: successor(RBNode x){
    if (x->right != nil){
        x = x->right;
        while (x->left != nil){
            x = x->left;
        }
    }else{
        while (x->parent != nil && x->isright()){
            x = x->parent;
        }
    }
    return x;
}

template <class T>
RBNode RBTree<T>:: predecessor(RBNode x){
    if (x->left != nil){
        x = x->left;
        while (x->right != nil){
            x = x->right;
        }
    }else{
        while (x->parent != nil && x->isleft()){
            x = x->parent;
        }
    }
    return x;
}

template<class T>
RBNode RBTree<T>:: get_lower_bound(RBNode x,T v){
    if (x == nil || x->val == v){
        return x;
    }else if (v > x->val){
        return get_lower_bound(x->right,v);
    }else{
        RBNode y = get_lower_bound(x->left,v);
        if (y != nil) return y;
    }
    return x;
}

template<class T>
RBNode RBTree<T>:: lower_bound(T x){
    return get_lower_bound(root,x);
}

template<class T>
RBNode RBTree<T>:: get_upper_bound(RBNode x,T v){
    if (x == nil){
        return x;
    }else if (v <= x->val){
        return get_upper_bound(x->left,v);
    }else{
        RBNode y = get_upper_bound(x->right,v);
        if (y != nil) return y;
    }
    return x;
}

template<class T>
RBNode RBTree<T>:: upper_bound(T x){
    return get_upper_bound(root,x);
}

template<class T>
inline RBNode RBTree<T>:: rfind(T v){
    // 直接调用 get_lower_bound 虽然方便，但是效率会略低
    // return get_lower_bound(root,v);
    RBNode x = root;
    while (x != nil){
        if (v > x->val) x = x->right;
        else if (v < x->val) x = x->left;
        else return x;
    }
    return x;
}

template<class T>
inline RBNode RBTree<T>:: findkth(int v){
    RBNode x = root;
    while (x != nil){
        if (1 + x->left->s > v){
            x = x->left;
        }else if (1 + x->left->s < v){
            v -= 1 + x->left->s;
            x = x->right;
        }else{
            return x;
        }
    }
    return x;
}

// 返回小于 v 的数的个数 + 1
template<class T>
inline int RBTree<T>:: get_rank(T v){
    RBNode x = root;
    int ret = 1;
    while (x != nil){
        if (v > x->val){
            ret += 1 + x->left->s;
            x = x->right;
        }else{
            x = x->left;
        }
    }
    return ret;
}

template<class T>
void RBTree<T>:: left_rotate(RBNode x){
    RBNode y = x->right;
    if (x == root){
        root = y;
    }else if (x->isleft()){
        x->parent->left = y;
    }else{
        x->parent->right = y;
    }
    y->parent = x->parent;
    x->right = y->left;
    if (y->left != nil) y->left->parent = x;
    y->left = x;
    x->parent = y;
    if (x != nil) x->update();
    if (y != nil) y->update();
}

template<class T>
void RBTree<T>:: right_rotate(RBNode x){
    RBNode y = x->left;
    if (x == root){
        root = y;
    }else if (x->isright()){
        x->parent->right = y;
    }else{
        x->parent->left = y;
    }
    y->parent = x->parent;
    x->left = y->right;
    if (y->right != nil) y->right->parent = x;
    y->right = x;
    x->parent = y;
    if (x != nil) x->update();
    if (y != nil) y->update();
}

template<class T>
void RBTree<T>:: insert(T v){
    RBNode x = root;
    RBNode y = nil;
    while (x != nil){
        y = x;
        if (v > x->val){
            x = x->right;
        }else{
            x = x->left;
        }
    }

    x = new Node<T>(y,nil,nil,1,v,RED);
    if (y == nil) root = x;
    else if (v > y->val){
        y->right = x;
    }else{
        y->left = x;
    }
    while (y != nil){
        y->update();
        y = y->parent;
    }
    FixDoubleRed(x);
}

template<class T>
void RBTree<T>:: FixDoubleRed(RBNode x){
    while (x != root && x->parent->color == RED){
        if (x->parent->isleft()){
            RBNode uncle = x->parent->parent->right;
            if (uncle->color == RED){
                x->parent->color = uncle->color = BLACK;
                x->parent->parent->color = RED;
                x = x->parent->parent; // 递归向上
            }else{
                if (x->isright()){ // 三角形变直线型
                    x = x->parent;
                    left_rotate(x);
                }
                x->parent->color = BLACK;
                x->parent->parent->color = RED;
                right_rotate(x->parent->parent);
                break;
            }
        }else{
            RBNode uncle = x->parent->parent->left;
            if (uncle->color == RED){
                x->parent->color = uncle->color = BLACK;
                x->parent->parent->color = RED;
                x = x->parent->parent; // 递归向上
            }else{
                if (x->isleft()){ // 三角形变直线型
                    x = x->parent;
                    right_rotate(x);
                }
                x->parent->color = BLACK;
                x->parent->parent->color = RED;
                left_rotate(x->parent->parent);
                break;  
            }
        }
    }
    root->color = BLACK;
}

template<class T>
void RBTree<T>:: remove(RBNode x){
    if (x == root && x->left == nil && x->right == nil){
        root = nil;
        delete x;
        return ;
    }

    RBNode y;
    if (x->left != nil && x->right != nil){
        y = successor(x);
        x->val = y->val;
        x = y;
    }
    // 现在 x 的两个儿子中至少有一个 nil
    if (x->color == RED){
        //不用管
    }else{
        if (x->right != nil){ // 则一定 x->left == nil, 且 x->right 为红色
            x->val = x->right->val;
            x = x->right;
        }else if (x->left != nil){
            x->val = x->left->val;
            x = x->left;
        }else{ // x->left == x->right == nil
            FixDoubleBlack(x);
        }
    }
    y = x->parent;
    if (x->isleft()) y->left = nil;
    else y->right = nil;
    while (y != nil){
        y->update();
        y = y->parent;
    }
    delete x;
}

template<class T>
void RBTree<T>:: FixDoubleBlack(RBNode x){
    while (x != root && x->color == BLACK){
        if (x->isleft()){
            RBNode brother = x->parent->right;
            if (brother->color == RED){
                x->parent->color = RED;
                brother->color = BLACK;
                left_rotate(x->parent);
                brother = x->parent->right; // 转换为 brother->color == BLACK
            }
            if (brother->left->color == BLACK && brother->right->color == BLACK){
                brother->color = RED;
                if (x->parent->color == RED){
                    x->parent->color = BLACK;
                    break;
                }else{
                    x = x->parent;
                }
            }else{
                if (brother->left->color == RED){
                    brother->left->color = BLACK;
                    brother->color = RED;
                    right_rotate(brother);
                    brother = x->parent->right;
                }
                brother->right->color = BLACK;
                brother->color = x->parent->color;
                x->parent->color = BLACK;
                left_rotate(x->parent);
                break;
            }
        }else{
            RBNode brother = x->parent->left;
            if (brother->color == RED){
                x->parent->color = RED;
                brother->color = BLACK;
                right_rotate(x->parent);
                brother = x->parent->left; // 转换为 brother->color == BLACK
            }
            if (brother->right->color == BLACK && brother->left->color == BLACK){
                brother->color = RED;
                if (x->parent->color == RED){
                    x->parent->color = BLACK;
                    break;
                }else{
                    x = x->parent;
                }
            }else{
                if (brother->right->color == RED){
                    brother->right->color = BLACK;
                    brother->color = RED;
                    left_rotate(brother);
                    brother = x->parent->left;
                }
                brother->left->color = BLACK;
                brother->color = x->parent->color;
                x->parent->color = BLACK;
                right_rotate(x->parent);
                break;
            }
        }
    }
    root->color = BLACK;
}

int main(){
    // luogu P6136
    
    RBTree<int> tree;
    int n,m,op,x,sum = 0,last = 0;  
    ios::sync_with_stdio(false);
    cin >> n >> m;
    for (int i=1;i<=n;++i){
        cin >> x;
        tree.insert(x);
    }
    for (int i=1;i<=m;++i){
        cin >> op >> x;
        x ^= last;
        //cout << (i+1) << ":  "; 
        //cout << "x = " << x << endl;
        if (op==1)       tree.insert(x);
        else if(op==2)   tree.remove(tree.rfind(x));
        else if(op==3)   last = tree.get_rank(x);
        else if (op==4)  last = tree.findkth(x)->val;
        else if (op==5)  last = tree.upper_bound(x)->val;
        else             last = tree.lower_bound(x+1)->val;
        //cout << "last = " << last << endl;
        if (op>2) sum ^= last;
        //tree.printStruct();
        
    }
    cout << sum;
    return 0;
}
