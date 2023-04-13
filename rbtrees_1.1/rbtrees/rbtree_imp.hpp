/**
 * @file rbtree_imp.hpp
 *
 * CopyRight F. J. Madrid-Cuevas <fjmadrid@uco.es>
 *
 * Sólo se permite el uso de este código en la docencia de las asignaturas sobre
 * Estructuras de Datos de la Universidad de Córdoba.
 *
 * Está prohibido su uso para cualquier otro objetivo.
 */
#pragma once
#include <algorithm>
#include <cmath>
#include <cassert>
#include <limits>

#include "rbtree.hpp"

/****
 * RBTNode class implementation.
 ****/

template <class T>
RBTNode<T>::RBTNode(T const &it,
                    RBTNode<T>::Color c,
                    RBTNode<T>::Ref p,
                    RBTNode<T>::Ref l,
                    RBTNode<T>::Ref r)
{
    // TODO
    this->item_ = it;
    this->color_ = c;
    this->parent_ = p;
    this->left_ = l;
    this->right_ = r;
    
    //
    assert(item() == it);
    assert(color() == c);
    assert(parent() == p);
    assert(left() == l);
    assert(right() == r);
}

template <class T>
typename RBTNode<T>::Ref RBTNode<T>::create(T const &it,
                                            RBTNode<T>::Color color,
                                            RBTNode<T>::Ref parent,
                                            RBTNode<T>::Ref left,
                                            RBTNode<T>::Ref right)
{
    auto ref = RBTNode<T>::Ref(new RBTNode<T>(it, color, parent, left, right));
    ref->this_ = ref;
    return ref;
}

template <class T>
RBTNode<T>::~RBTNode()
{
    // TODO
    // Hint: Think, is it needed to do anything?.

    //
}

template <class T>
T RBTNode<T>::item() const
{
    T value;
    // TODO
    value = this->item_;
     //
    return value;
}

template <class T>
typename RBTNode<T>::Ref RBTNode<T>::parent() const
{
    RBTNode<T>::Ref node;
    // TODO
    node = this->parent_;
    //
    return node;
}

template <class T>
typename RBTNode<T>::Ref RBTNode<T>::left() const
{
    RBTNode<T>::Ref node;
    // TODO
    node = this->left_;
    //
    return node;
}

template <class T>
typename RBTNode<T>::Ref RBTNode<T>::right() const
{
    RBTNode<T>::Ref node;
    // TODO
    node = this->right_;
    //
    return node;
}

template <class T>
typename RBTNode<T>::Ref RBTNode<T>::child(int idx) const
{
    assert(idx == 0 || idx == 1);
    RBTNode<T>::Ref node;
    // TODO
    return (idx == 0) ? right() : left();
    //
    assert(idx == 0 || node == right());
    assert(idx == 1 || node == left());
    return node;
}

template <class T>
typename RBTNode<T>::Color RBTNode<T>::color() const
{
    RBTNode<T>::Color color = BLACK;
    // TODO
    color = this->color_;
    //
    return color;
}

template <class T>
const typename RBTNode<T>::Ref &RBTNode<T>::This() const
{
    return this_;
}

template <class T>
typename RBTNode<T>::Ref RBTNode<T>::This()
{
    return this_;
}

template <class T>
void RBTNode<T>::set_item(const T &new_it)
{
    // TODO
    this->item_ = new_it;
    //
    assert(item() == new_it);
}

template <class T>
void RBTNode<T>::set_color(RBTNode<T>::Color new_color)
{
    // TODO
    this->color_ = new_color;
    //
    assert(color() == new_color);
}

template <class T>
void RBTNode<T>::set_parent(RBTNode<T>::Ref new_parent)
{
    // TODO
    this->parent_ = new_parent;
    //
    assert(parent() == new_parent);
}

template <class T>
void RBTNode<T>::set_left(Ref new_child)
{
    // TODO
    // Update the right child link
    left_ = new_child;

    // Update the parent link of the new right child
    if (new_child) {
        new_child->parent_ = This();
    }
    assert(left() == new_child);
    assert(!new_child || new_child->parent() == This());
}

template <class T>
void RBTNode<T>::set_right(RBTNode<T>::Ref new_child)
{
    // TODO
    // Remember to update the parent link of new_child too.
    // Update the right child link
    right_ = new_child;

    // Update the parent link of the new right child
    if (new_child) {
        new_child->parent_ = This();
    }
    //
    assert(right() == new_child);
    assert(!new_child || new_child->parent() == This());
}

template <class T>
void RBTNode<T>::set_child(int idx, RBTNode<T>::Ref new_child)
{
    assert(idx == 0 || idx == 1);
    // TODO
    // Remember to update the parent link of new_child too.
    if (idx == 0) {
        left_ = new_child;
    } else {
        right_ = new_child;
    }
    // Update parent link of new_child
    if (new_child) {
        new_child->parent_ = This();
    }
    //
    assert(idx == 0 || new_child == right());
    assert(idx == 1 || new_child == left());
    assert(!new_child || new_child->parent() == This());
}

/***
 * RBTree class implementation.
 ***/

template <class T>
RBTree<T>::RBTree()
{
    // TODO
    this->root_ = nullptr;
    this->parent_ = nullptr;
    this->current_= nullptr;
    //
    assert(is_a_binary_search_tree());
    assert(is_a_rbtree());
    assert(is_empty());
}

template <class T>
RBTree<T>::RBTree(T const &item)
{
    // TODO
    this->root_ = RBTNode<T>::create(item, RBTNode<T>::BLACK, nullptr, nullptr, nullptr);
    //
    assert(is_a_binary_search_tree());
    assert(is_a_rbtree());
    assert(!is_empty());
}

template <class T>
const typename RBTree<T>::Ref &RBTree<T>::This() const
{
    return this_;
}

template <class T>
typename RBTree<T>::Ref RBTree<T>::This()
{
    return this_;
}

template <class T>
typename RBTree<T>::Ref RBTree<T>::create()
{
    auto ref = Ref(new RBTree<T>());
    ref->this_ = ref;
    return ref;
}

template <class T>
typename RBTree<T>::Ref RBTree<T>::create(T const &item)
{
    auto ref = Ref(new RBTree<T>(item));
    ref->this_ = ref;
    return ref;
}

template <class T>
typename RBTree<T>::Ref RBTree<T>::create(std::istream &in) noexcept(false)
{
    auto tree = RBTree<T>::create();
    std::string token;
    // TODO
    // Hint: use the recursive definition of a tree to unfold.
    if (token == "[") //En el caso en el que el stream no empiece por [, la única opción es que comience por []
    {
        T data; //En este item almacenamos el item que leamos del stream y que meteremos en el árbol
        in >> data;
        tree = RBTree<T>::create(data);

        char color;
        in >> color;
        if (color == 'R')
            tree->current_->set_color(RBTNode<T>::RED);
        else if (color == 'B')
            tree->current_->set_color(RBTNode<T>::BLACK);
        else
            throw std::runtime_error("Wrong input format");

        auto l_subtree = RBTree<T>::create(in);
        auto r_subtree = RBTree<T>::create(in);
        tree->set_left(l_subtree);
        tree->set_right(r_subtree);

        in >> token;
    }
    else if (token != "[]") //Sí el stream representa un árbol vacío y el formato no es [], nos da error de formato
        throw std::runtime_error("Wrong input format");
    //
    tree->this_ = tree;
    if (!tree->is_a_binary_search_tree())
        throw std::runtime_error("It is not a binary search tree");
    else if (!tree->is_a_rbtree())
        throw std::runtime_error("It is not an rbtree");
    return tree;
}

template <class T>
RBTree<T>::~RBTree()
{
    // TODO
    // Hint: Think, is it needed to do anything?.

    //
}
#ifdef __ONLY_BSTREE__
/**
 * @brief Create an perfectly balanced BSTree by inserting the median of
 *        an ordered sub sequence data[begin, end).
 * @param data is the ordered sequence of values.
 * @param begin,
 * @param end specify a [begin, end) indexing interval of array data to use.
 * @pre 0 <= begin <= end <=data.size()
 * @pre begin==end || data[begin]<data[end];
 */
template <class T>
void create_inserting_median(std::vector<T> const &data,
                             size_t begin,
                             size_t end,
                             typename RBTree<T>::Ref &tree)
{
    assert(begin <= end);
    assert(end <= data.size());
    assert(begin == end || data[begin] <= data[end - 1]);

    // TODO
    // Hint: if (end==begin) none thing must be done (it is an empty sub array)
    //  else, insert the median in the tree and (recursively) process
    //  the two sub sequences [begin, median_idx) and [median_idx+1, end)
    if ((end - begin) >= 1) //Entramos siempre que la resta mayor o igual a 1
    {
        tree->insert(data[begin + ((end - begin) / 2)]);
        //Llamamos recursivamente a la funcion
        create_inserting_median(data, begin, begin + ((end - begin) / 2), tree);
        //Llamamos recursivamente a la funcion
        create_inserting_median(data, (begin + ((end - begin) / 2)) + 1, end, tree);
    }
    //
}

template <class T>
typename RBTree<T>::Ref
create_perfectly_balanced_bstree(std::vector<T> &data)
{
    typename RBTree<T>::Ref tree = RBTree<T>::create();
    // TODO
    // Remember: the empty tree is perfectly balanced.
    // Remember: first, an ordered sequence (using < order) of values is needed.
    // Then you should use the above create_inserting_median function
    // on a empty tree to recursively create the perfectly balanced bstree.
    std::sort(data.begin(), data.end());
    create_inserting_median(data, 0, data.size(), tree);
    //
    assert(tree != nullptr);
    return tree;
}
#endif //__ONLY_BSTREE__

template <class T>
bool RBTree<T>::is_empty() const
{
    bool empty = false;
    // TODO
    if(this->root_ == nullptr){
        empty = true;
    }
    //
    return empty;
}

template <class T>
T RBTree<T>::item() const
{
    assert(!is_empty());
    T value;
    // TODO
    value = this->root_->item();
    //
    return value;
}

template <class T>
std::ostream &RBTree<T>::fold(std::ostream &out) const
{
    // TODO
    // Hint: use the recursive definition of a tree to fold it.
    if (this->root_ == nullptr) //Caso en el que el árbol está vacío. El output será []
         out << "[]";

     else //Caso en el que el árbol no está vacío
     {
         out << "[ ";
         out << this->item();
         out << " ";
         this->left()->fold(out);
         out << " ";
         this->right()->fold(out);
         out << " ]";
     }

    //
    return out;
}

template <class T>
bool RBTree<T>::current_exists() const
{
    bool exists = false;
    // TODO
    exists = this->current_ != nullptr;
    //
    return exists;
}

template <class T>
T RBTree<T>::current() const
{
    assert(current_exists());
    T value;
    // TODO
    value = this->current_->item();
    //
    return value;
}

template <class T>
int RBTree<T>::current_level() const
{
    assert(current_exists());
    int level = 0;
    // TODO
    // Hint: follow the chain of parents' links.
    auto node_aux = this->root_; //Creamos un nodo que tendra el valor de root_ y nos servira de ayuda
    while (this->current_ != node_aux)
    {
        if (node_aux->item() > this->current_->item()) //Entramos si el item del nodo es mayor que el item de current
            node_aux = node_aux->left();

        else if (this->current_->item() > node_aux->item()) //Entramos si el item del nodo es menor que el item de current
            node_aux = node_aux->right();

        level++;
    }
    //
    return level;
}

template <class T>
typename RBTree<T>::Ref RBTree<T>::left() const
{
    assert(!is_empty());
    typename RBTree<T>::Ref subtree = nullptr;
    // TODO
    // Hint: use the protected method to create a tree given the root node.
    if (root_->left()) {  // If left child exists, create a new subtree from it
        subtree = RBTree<T>::create(root_->left());
    }    if (this->root_->left() != nullptr)
    {
        subtree->create_root(root_->left()->item());
        subtree->root_->set_parent(nullptr);
        subtree->root_->set_left(root_->left()->left());
        subtree->root_->set_right(root_->left()->right());
    }
    //
    return subtree;
}

template <class T>
typename RBTree<T>::Ref RBTree<T>::right() const
{
    assert(!is_empty());
    typename RBTree<T>::Ref subtree = nullptr;
    // TODO
    // Hint: use the protected method to create a tree given the root node.
    if (this->root_->right() != nullptr)
    {
        subtree->create_root(root_->right()->item());
        subtree->root_->set_parent(nullptr);
        subtree->root_->set_left(root_->right()->left());
        subtree->root_->set_right(root_->right()->right());
    }
    //
    return subtree;
}

template <class T>
size_t compute_size(typename RBTree<T>::Ref const &tree)
{
    assert(tree != nullptr);
    int s = 0;
    // TODO
    // Hint: use the recursive implementation.
    if (tree == nullptr)
        return 0;

    size_t left_size = compute_size<T>(tree->left());
    size_t right_size = compute_size<T>(tree->right());
    s = left_size + right_size + 1;
    //
    return s;
}

template <class T>
int compute_height(typename RBTree<T>::Ref const &tree)
{
    assert(tree != nullptr);
    int h = -1;
    // TODO
    // Hint: use the recursive implementation.
    if (!tree->is_empty()) {
        int left_height = compute_height<T>(tree->left());
        int right_height = compute_height<T>(tree->right());
        h = std::max(left_height, right_height) + 1;
    }
    //
    return h;
}

template <class T>
bool RBTree<T>::has(const T &k) const
{
#ifndef NDEBUG
    bool old_current_exists = current_exists();
    T old_current;
    if (old_current_exists)
        old_current = current();
#endif

    bool found = true;

    // TODO
    // Hint: you can reuse the search method for this but in this case you will
    //       need to use "const_cast" to remove constness of "this" and
    //       save/restore the old state of current before returning.
    // Search for the key in the tree
    auto node_aux = this->root_;
    int i = 1;
    while (i>0) //Bucle infinito utilizamos los break para salir de el en caso de que hayamos terminado
    {
        if (node_aux->item() > k) //k esta por la izquierda
        {
            if (node_aux->left() != nullptr) //SI el nodo tiene izquierdo entramos
                node_aux = node_aux->left();
            else
            {
                found = false;
                break;
            }
        }
        else if (node_aux->item() < k) //k esta por la derecha
        {
            if (node_aux->right() != nullptr) //Si el nodo tiene derecho entramos
                node_aux = node_aux->right();
            else
            {
                found = false;
                break;
            }
        }
        else
            break;
    }
    //
#ifndef NDEBUG
    assert(!old_current_exists || old_current == current());
#endif
    return found;
}

/**
 * @brief infix process of a node.
 * The Processor must allow to be used as a function with a parameter  (the
 * item to be processed) and returning true if the process must continue or
 * false if not.
 * @param node is the node to be processed.
 * @param p is the Processor.
 * @return true if all the tree was in-fix processed.
 */
template <class T, class Processor>
bool infix_process(typename RBTNode<T>::Ref node, Processor &p)
{
    bool retVal = true;
    // TODO
    // Remember: if node is nullptr return true.
    if (!node) {
        return retVal;
    }

    if (node->left && !infix_process<T, Processor>(node->left, p)) {
        retVal=false;
    }

    if (!p(node->item)) {
        retVal=false;
    }

    if (node->right && !infix_process<T, Processor>(node->right, p)) {
        retVal=false;
    }

    //
    return retVal;
}

template <class T>
bool RBTree<T>::is_a_binary_search_tree() const
{
    bool is_bst = true;
    // TODO
    // Remember: a empty tree is a binary search tree.
    //
    // Remember: for a non empty binary search tree, the in-fix traversal from
    // the root node must follow an ordered sequence of items.
    //
    // Remember: use a lambda function with signature '(T v) -> bool' to
    //  implement the Processor.
    //
    // Lambda function to check if a value is greater than or equal to a given value

    //
    return is_bst;
}

template <class T>
bool RBTree<T>::meet_red_invariant() const
{
#ifdef __ONLY_BSTREE__
    return true;
#else
    bool is_met = true;
    // TODO
    // Remember: A red node must not have a red child.
    // Remember: An empty tree meets the invariant.

    //
    return is_met;
#endif
}

template <class T>
bool RBTree<T>::meet_black_invariant() const
{
#ifdef __ONLY_BSTREE__
    return true;
#else
    bool is_met = true;
    // TODO
    // Remember: The number of black nodes for each path from root to any leaf must be equal.
    // Hint: use a lambda function to travel the tree.
    // Check red nodes for having red children

    //
    return is_met;
#endif
}

template <class T>
bool RBTree<T>::is_a_rbtree() const
{
#ifdef __ONLY_BSTREE__
    return true;
#else
    assert(meet_red_invariant());
    assert(meet_black_invariant());
    return true;
#endif
}

template <class T>
void RBTree<T>::create_root(T const &v)
{
    assert(is_empty());
    // TODO
    this->root_ = RBTNode<T>::create(v, RBTNode<T>::BLACK ,nullptr, nullptr,nullptr);
    //
    assert(is_a_binary_search_tree());
    assert(is_a_rbtree());
    assert(!is_empty());
    assert(item() == v);
}

template <class T>
bool RBTree<T>::search(T const &k)
{
    bool found = false;
    // TODO
    this->parent_ = nullptr;
    this->current_ = this->root_;
    while (this->current_ != nullptr) //Mientras current no sea nulo y no encontremos el valor entramos en el bucle
    {
        if (this->current_->item() == k) //Entramos en el bucle si encontramos el nodo
        {
            found = true;
            break;
        } //Si lo encontramos saldremos del bucle al acabar
        else
        {
            this->parent_ = this->current_; //EL padre sera nuestro actual current, ya que current pasara a ser o su hijo izq o dcho
            if (this->current_->item() < k) //Si el tenemos es menor que el que buscamos, tenemos que ir a la derecha
                this->current_ = this->current_->right();
            else //Si el que tenemos es mayor que el que buscamos, tenemos que ir a la izquierda
                this->current_ = this->current_->left();
        }
    }
    //
    assert(!found || current() == k);
    assert(found || !current_exists());
    return found;
}

template <class T>
void RBTree<T>::insert(T const &k)
{
    // //Check invariants.
    assert(is_a_binary_search_tree());
    assert(is_a_rbtree());

    if (!search(k))
    {
        // TODO
        // Remember: a new node is always RED.
        // Create a new node as a red node.
        auto new_node = RBTNode<T>::create(k, RBTNode<T>::RED ,nullptr, nullptr,nullptr);

        if (this->is_empty())
        {
            // If the tree is empty, make the new node as the root node.
            this->root_ = std::move(new_node);
        }
        else
        {
            // Otherwise, find the correct position to insert the new node.
            auto current_node = this->root_;
            while (true)
            {
                if (k < current_node->item())
                {
                    if (current_node->left() != nullptr)
                    {
                        current_node = current_node->left();
                    }
                    else
                    {
                        current_node->set_left(std::move(new_node));
                        break;
                    }
                }
                else
                {
                    if (current_node->right() != nullptr)
                    {
                        current_node = current_node->right();
                    }
                    else
                    {
                        current_node->set_right(std::move(new_node));
                        break;
                    }
                }
            }
        }

        //

        assert(check_parent_chains());
        make_red_black_after_inserting();
        assert(check_parent_chains());
        assert(check_min_max_branch_length());
    }
    //

    // Check invariants.
    assert(is_a_binary_search_tree());
    assert(is_a_rbtree());

    // check postconditions.
    assert(current_exists());
    assert(current() == k);
}

template <class T>
void RBTree<T>::remove()
{
    // check preconditions.
    assert(current_exists());

#ifndef NDEBUG
    // the invariants only must be checked for the first recursive call.
    // We use a static variable to count the recursion levels.
    // see section "Static variables in a Function" in
    // Ref https://www.geeksforgeeks.org/static-keyword-cpp/ for more info.
    static int recursion_count = 0;
    recursion_count++;
    if (recursion_count == 1)
    {
        // Check invariants.
        assert(is_a_binary_search_tree());
        assert(is_a_rbtree());
    }
#endif // NDEBUG

    bool replace_with_subtree = true;
    typename RBTNode<T>::Ref subtree;

    // TODO
    //  Check which of BSTree cases 0,1,2,3 we have (see theoretical class slides).
    if (!this->current_->right() && !this->current_->left()) //Si no tenemos ni izq ni dcho, entramos
            subtree = nullptr;

        else if (!this->current_->left()) //Entramos si no tiene izq
            subtree = this->current_->right();

        else if (!this->current_->right()) //Entramos si no tiene dcho
            subtree = this->current_->left();

        else //Entramos si tenemos izquierdo y derecho
            replace_with_subtree = false;

    //

    if (replace_with_subtree)
    {
        // TODO
        // Manage cases 0,1,2
        if (this->parent_ == nullptr) //Entramos si el padre es nulo
            this->root_ = subtree;

        else if (this->current_ == this->parent_->left()) //Entramos si el actual es igual al izq del padre
            this->parent_->set_left(subtree);

        else if (this->current_ == this->parent_->right()) //Entramos si el actual es igual al dcho del padre
            this->parent_->set_right(subtree);

        this->current_ = nullptr;
        //
        assert(check_parent_chains());
        make_red_black_after_removing(subtree);
        assert(check_parent_chains());
        assert(check_min_max_branch_length());

        set_current_node(nullptr);
    }
    else
    {
        // TODO
        // Manage case 3.
        auto x = this->current_;
        find_inorder_successor();
        x->set_item(this->current_->item());
        remove();
        //
    }

#ifndef NDEBUG
    // We come back so the recursion count must be decreased.
    recursion_count--;
    assert(recursion_count >= 0);
    if (recursion_count == 0)
    {
        // Only check for the last return.
        // Check invariants.
        assert(is_a_binary_search_tree());
        assert(is_a_rbtree());

        // Check postconditions.
        assert(!current_exists());
    }
#endif
}

template <class T>
RBTree<T>::RBTree(typename RBTNode<T>::Ref root_node)
{
    // TODO
    // Set the root node of the tree
    root_ = root_node;

    // If there is a root node, set its parent to nullptr and ensure it is black
    if (root_ != nullptr) {
        root_->set_parent(nullptr);
        root_->set_color(RBTNode<T>::BLACK);
    }

    // Set the current node to nullptr
    set_current_node(nullptr);
    //
    assert(!current_exists());
}

template <class T>
typename RBTree<T>::Ref RBTree<T>::create(typename RBTNode<T>::Ref root)
{
    // Note: we can not use std:make_shared here because the
    //  constructor with a root node is private.
    Ref ref = Ref(new RBTree<T>(root));
    ref->this_ = ref;
    return ref;
}

template <class T>
void RBTree<T>::set_left(Ref subtree)
{
    assert(!is_empty());

    // TODO
    if (subtree->is_empty()) //Si le pasamos un arbol vacio, seteamos nullptr
        this->root_->set_left(nullptr);

    else
        this->root_->set_left(RBTNode<T>::create(subtree->root_->item(),subtree->root_->color(), this->root_, subtree->root_->left(), subtree->root_->right()));
    //

    assert(subtree->is_empty() || left()->item() == subtree->item());
    assert(!subtree->is_empty() || left()->is_empty());
}

template <class T>
void RBTree<T>::set_right(Ref subtree)
{
    assert(!is_empty());

    // TODO
    if (subtree->is_empty()) //Si le pasamos un arbol vacio, seteamos nullptr
        this->root_->set_right(nullptr);

    else
        this->root_->set_right(RBTNode<T>::create(subtree->root_->item(),subtree->root_->color(), this->root_, subtree->root_->left(), subtree->root_->right()));
    //
    assert(subtree->is_empty() || right()->item() == subtree->item());
    assert(!subtree->is_empty() || right()->is_empty());
}

template <class T>
typename RBTNode<T>::Ref RBTree<T>::current_node() const
{
    typename RBTNode<T>::Ref node = nullptr;
    // TODO
    if (current_ != nullptr) {
        node = current_;
    }
    //
    return node;
}

template <class T>
void RBTree<T>::set_current_node(typename RBTNode<T>::Ref new_curr)
{
    // TODO
    current_ = new_curr;
    //
    assert(new_curr == current_node());
}

template <class T>
typename RBTNode<T>::Ref RBTree<T>::root_node() const
{
    typename RBTNode<T>::Ref node = nullptr;
    // TODO
    if (!is_empty()) {
        node = current_node();
        while (node->parent() != nullptr) {
            node = node->parent();
        }
    }
    //
    return node;
}

template <class T>
void RBTree<T>::set_root_node(typename RBTNode<T>::Ref new_root)
{
    // TODO
    this->root_ = new_root;
    //
    assert(new_root == root_node());
}

template <class T>
typename RBTNode<T>::Ref RBTree<T>::parent_node() const
{
    typename RBTNode<T>::Ref node = nullptr;
    // TODO
    if (current_exists()) {
        node = current_node()->parent();
    }
    //
    return node;
}

template <class T>
void RBTree<T>::find_inorder_successor()
{
    assert(current_exists());
    assert(is_a_binary_search_tree());
#ifndef NDEBUG
    T old_curr = current();
#endif
    // TODO
    if (current_node()->right() != nullptr) {
        // Case 1: The right subtree is not empty.
        set_current_node(current_node()->right());
        while (current_node()->left() != nullptr) {
            set_current_node(current_node()->left());
        }
    } else {
        // Case 2: The right subtree is empty.
        typename RBTNode<T>::Ref parent = parent_node();
        while (parent != nullptr && current_node() == parent->right()) {
            set_current_node(parent);
            parent = parent_node();
        }
        set_current_node(parent);
    }
    //
    assert(current_exists() && current_node()->left() == nullptr);
#ifndef NDEBUG
    assert(current() > old_curr);
#endif
}

template <class T>
typename RBTNode<T>::Ref
RBTree<T>::rotate(typename RBTNode<T>::Ref P,
                  int dir)
{
    assert(P != nullptr);
    assert(dir == 0 || dir == 1);
    assert(P->child(1 - dir) != nullptr);
    typename RBTNode<T>::Ref N = P->child(1 - dir); // the child to promote.
#ifdef __DEBUG__
    if (__DEBUG__ > 1)
        std::clog << "Rotating to " << (dir == 0 ? "left" : "right") << " on key " << P->item() << std::endl;
#endif
    // TODO
    // Remember update links to parents.
    // Hint: you can see wikipedia: https://en.wikipedia.org/wiki/Tree_rotation
    typename RBTNode<T>::Ref CN = N->child(dir); // the child to promote.
    typename RBTNode<T>::Ref G; // Grandparent.
    P->set_child(1-dir,CN);
    N->set_child(dir,P);
    if(G == nullptr){
        G->set_child(dir,N);
    }else
        this->set_root_node(N);
    //
    return N;
}

template <class T>
void RBTree<T>::make_red_black_after_inserting()
{
#ifdef __ONLY_BSTREE__
    return;
#else
    // TODO
    // Hint: @see https://www.geeksforgeeks.org/insertion-in-red-black-tree/

    typename RBTNode<T>::Ref N = current_node();
    typename RBTNode<T>::Ref P = parent_node();
    typename RBTNode<T>::Ref G; // Grandparent.
    typename RBTNode<T>::Ref U; // Uncle.
    int g_p_dir; //Granddad to parent direction.
    int p_n_dir; //Parent to N direction.

    if (P == nullptr)
    {
        //TODO: Case 0:
        N->set_color(RBTNode<T>::BLACK);
        set_root_node(N);
        return;
        //
    }

    while (P != nullptr)
    {
        if (P->color() == RBTNode<T>::BLACK)
        {
            //TODO Case 1: reqs 3 and 4 are met.
            return;
            //
        }

        // From here P is red.

        //TODO: update G, g_p_dir, U.
        G = P->parent();
        g_p_dir = (G->left() == P) ? 0 : 1;
        p_n_dir = (P->left() == N) ? 0 : 1;
        U = (g_p_dir == 0) ? G->right() : G->left();
        //

        if (U != nullptr && U->color() == RBTNode<T>::RED)
        {
            //TODO Case 2:
            P->set_color(RBTNode<T>::BLACK);
                       U->set_color(RBTNode<T>::BLACK);
                       G->set_color(RBTNode<T>::RED);
                       // Set N to G, and update P, G, and U for next iteration.
                       N = G;
                       P = N->parent();
                       G = P->parent();
                       g_p_dir = (G->left() == P) ? 0 : 1;
                       U = (g_p_dir == 0) ? G->right() : G->left();
            //
        }
        else
        {
            //Case 3:

            //TODO: update p_n_dir

            //
            if (g_p_dir != p_n_dir)
            {
                // TODO: cases 3c (LR) 3d (RL)
                if (p_n_dir == 0)
                {
                    rotate(P,g_p_dir);
                    // Update necessary nodes and directions
                    N = N->left();
                    g_p_dir = 0;
                }
                else
                {
                    rotate(P,g_p_dir);
                    // Update necessary nodes and directions
                    N = N->right();
                    g_p_dir = 1;
                }

                //
            }

            // TODO: cases 3a (LL) 3b (RR)
            if (g_p_dir == 0)
            {
                rotate(G,g_p_dir);
            }
            // Case 3b: Right-Right case
            else
            {
                rotate(G,g_p_dir);
            }
            //
            P->set_color(RBTNode<T>::BLACK);
            G->set_color(RBTNode<T>::RED);
            return;
        }
        P = N->parent();
    }
#endif
}
template <class T>
void RBTree<T>::make_red_black_after_removing(typename RBTNode<T>::Ref V)
{

#ifdef __ONLY_BSTREE__
    return;
#else
    typename RBTNode<T>::Ref N = current_node(); // The removed node.
    typename RBTNode<T>::Ref P = parent_node();
    assert(!P || (P->child(0)==V || P->child(1)==V)); //V replaced to N

    if (P == nullptr)
    {
        // Case 1: empty tree.
        N->set_color(RBTNode<T>::BLACK);
        //
        return;
    }
    assert(P->child(0) == V || P->child(1) == V);
    int p_n_dir = (P->child(0)==V) ? 0 : 1; // Direction from P to N (0 Left, 1 Right)    

    if (N->color() == RBTNode<T>::RED || (V != nullptr && V->color() == RBTNode<T>::RED))
    {
        // TODO: case 2 (N or V is red)
        N->set_color(RBTNode<T>::BLACK);
        //
        return;
    }

    
    typename RBTNode<T>::Ref S; // Sibling of N
    typename RBTNode<T>::Ref C; // Close nephew of N
    typename RBTNode<T>::Ref D; // Distant nephew of N

    // Case 3: N and V are black.
    while (P != nullptr)
    {
        // TODO update S, D, C according to p_n_dir
        // Remember: S could be void.
        S = P->child(1 - p_n_dir);
        C = S->child(p_n_dir);
        D = S->child(1 - p_n_dir);
        //

        if (S == nullptr || S->color() == RBTNode<T>::BLACK)
        {
            //Case 3.1: N, V and S are black

            if ((C && C->color() == RBTNode<T>::RED) ||
                (D && D->color() == RBTNode<T>::RED))
            {
                // Case 3.1a (at least one nephew is red)
                if (D == nullptr || D->color() == RBTNode<T>::BLACK)
                {
                    // TODO: case 3.1a when only C is red (RL, LR)
                    // Remember update new D and S.
                    if (C->color() == RBTNode<T>::RED) // RL
                                        rotate(S,1);
                                     rotate(P,0);

                                     if (N != nullptr)
                                         N->set_color(RBTNode<T>::BLACK);
                                     S->set_color(RBTNode<T>::RED);
                                     return;
                    //
                }
                // TODO: case 3.1a when D is Red (RR, LL)
                if (N != nullptr)
                    N->set_color(RBTNode<T>::BLACK);
                S->set_color(RBTNode<T>::RED);
                D->set_color(RBTNode<T>::BLACK);
                //
                return;
            }
            else
            {
                // Case 3.1b C y D are black.
                if (P->color() == RBTNode<T>::RED)
                {
                    // TODO: Case 3.1b (parent is RED)
                    P->set_color(RBTNode<T>::BLACK);
                    S->set_color(RBTNode<T>::RED);
                    //
                    return;
                }
                else
                {
                    // TODO: Case 3.1b (parent is black)
                    // Remember: we must go up one level, so
                    // update N, P, and p_n_dir according (if new P<>Void).
                    S->set_color(RBTNode<T>::RED);
                    N = P;
                    P = N->parent();
                    if (P)
                        p_n_dir = (P->child(0) == N) ? 0 : 1;
                    //
                }
                //
            }
        }
        else
        {
            // TODO: Case 3.2: N,V are black, S is red
            rotate(P,0);
            P->set_color(RBTNode<T>::RED);
            S->set_color(RBTNode<T>::BLACK);
            P->child(p_n_dir)->set_color(RBTNode<T>::RED);
            P = S;
            S = P->child(1 - p_n_dir);
            C = S->child(p_n_dir);
            D = S->child(1 - p_n_dir);
            //
        }
    }
    return;
#endif
}

template <class T>
bool RBTree<T>::check_parent_chains()
{
    if (!is_empty())
    {
        std::function<void(typename RBTNode<T>::Ref, std::vector<T>)> go_down;
        go_down = [&go_down](typename RBTNode<T>::Ref node, std::vector<T> branch) -> void
        {
            if (node->left() != nullptr || node->right() != nullptr)
            {
                branch.push_back(node->item());
                if (node->left())
                    // go down by the left
                    go_down(node->left(), branch);
                if (node->right())
                    // go down by the right
                    go_down(node->right(), branch);
            }
            else
            {
                // The node is a leaf node, so check the branch
                // to the tree root node.
                typename RBTNode<T>::Ref parent = node->parent();
                int idx = static_cast<int>(branch.size()) - 1;
                while (parent && idx >= 0)
                {
                    assert(parent->item() == branch[idx]);
                    --idx;
                    parent = parent->parent();
                }
                assert(idx == -1 && parent == nullptr);
            }
        };
        std::vector<T> branch;
        go_down(root_node(), branch);
    }
    return true;
}

template <class T>
std::tuple<size_t, size_t>
compute_min_max_branch_length(typename RBTree<T>::Ref tree)
{
    assert(tree != nullptr);
    size_t min_path_l = 0;
    size_t max_path_l = 0;
    // TODO
    //  Hint: you can use a lambda function with prototype:
    //  std::function<void(typename RBTree<T>::Ref subt, size_t & min, size_t & max, size_t curr_l)>
    //  to recursive go down through the tree.
    //  Each new recursion increases the level in the tree. When a subtree is empty
    //  you have achieved a leaf of a branch and the current level is the
    //  branch length. Update the min/max values according before returning from
    //  this recursion.
    //  See: https://stackoverflow.com/questions/2067988/recursive-lambda-functions-in-c11
    //  to study a similar case.

    std::function<void(typename RBTree<T>::Ref, size_t &, size_t &, size_t)> go_down;
    go_down = [&go_down](typename RBTree<T>::Ref subtree, size_t &min_path_l,
                         size_t &max_path_l, size_t current_level) -> void
    {
        assert(!subtree->is_empty());
        auto left = subtree->left();
        auto right = subtree->right();
        if (left->is_empty() && right->is_empty())
        {
            min_path_l = std::min(min_path_l, current_level);
            max_path_l = std::max(max_path_l, current_level);
        }
        else
        {
            if (!left->is_empty())
                go_down(left, min_path_l, max_path_l, current_level + 1);
            if (!right->is_empty())
                go_down(right, min_path_l, max_path_l, current_level + 1);
        }
    };
    if (!tree->is_empty())
    {
        min_path_l = std::numeric_limits<size_t>::max();
        max_path_l = 0;
        go_down(tree, min_path_l, max_path_l, 0);
    }
    //

    return std::make_tuple(min_path_l, max_path_l);
}

template <class T>
bool RBTree<T>::check_min_max_branch_length() const
{
#ifdef __ONLY_BSTREE__
    return true;
#else
    size_t min_path_l, max_path_l;
    std::tie(min_path_l, max_path_l) = compute_min_max_branch_length<T>(This());
    return (max_path_l + 1) <= (2 * (min_path_l + 1));
#endif
}
