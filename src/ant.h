#ifndef ANT_H
#define ANT_H

#include <cstdint>
#include <vector>
#include <cassert>
#include <numeric>
#include <algorithm>

/*
Ant class for use in the ACS algorithm.

A single ant is always located at some node and represents a partial
(possibly complete) solution to a problem
*/
class Ant {
public:

    Ant(uint32_t nodes_count) :
        nodes_count_(nodes_count),
        value_(std::numeric_limits<double>::max()),
        unvisited_(nodes_count),
        index_in_unvisited_(nodes_count),
        is_visited_(nodes_count_, false) {

        std::iota(unvisited_.begin(), unvisited_.end(), 0);
        std::iota(index_in_unvisited_.begin(), index_in_unvisited_.end(), 0);

        // The last one is a guard - if anyone wants to check if node with index
        // equal to nodes_count is visited it will return true
        is_visited_.push_back(true);
    }


    Ant(const Ant &other) :
        nodes_count_(other.nodes_count_),
        value_(other.value_),
        position_(other.position_),
        visited_(other.visited_),
        unvisited_(other.unvisited_),
        index_in_unvisited_(other.index_in_unvisited_),
        is_visited_(other.is_visited_)
    {}


    virtual void reset() {
        *this = Ant(nodes_count_);
    }


    virtual ~Ant() {}


    bool is_visited(uint32_t node) const { return is_visited_[node]; }


    virtual bool is_available(uint32_t node) const { return !is_visited(node); }


    const std::vector<uint32_t> &get_visited() const { return visited_; }


    double get_value() const noexcept { return value_; }


    void set_value(double value) { value_ = value; }


    uint32_t get_position() const noexcept { return position_; }

    /*
     * Moves the ant to to a new node.
     * The method has O(1) complexity.
     */
    virtual void move_to(uint32_t node) {
        assert( visited_.size() < nodes_count_ );
        assert( !is_visited(node) );

        position_ = node;
        visited_.push_back(node);
        is_visited_[node] = true;

        remove_from_unvisited(node);
    }


    std::vector<uint32_t> &get_unvisited() {
        unvisited_.resize(nodes_count_ - visited_.size());
        return unvisited_;
    }


    bool has_complete_solution() const {
        return std::find(is_visited_.begin(), is_visited_.end(), false) == is_visited_.end();
    }


    void set_visited(std::vector<uint32_t> visited) {
        // We only support replacing the complete solution
        // otherwise is_visited_ and other members should be adjusted to the
        // modified contents
        assert(has_complete_solution());

        visited_ = visited;
    }

protected:

    void remove_from_unvisited(uint32_t node) {
        auto tail_index = nodes_count_ - visited_.size();
        auto node_index = index_in_unvisited_[node];
        if (node_index < tail_index) {
            auto tail = unvisited_[tail_index];
            unvisited_[node_index] = tail;
            index_in_unvisited_[tail] = node_index;
        }
    }

    uint32_t nodes_count_;
    double value_;
    uint32_t position_ = 0;
    std::vector<uint32_t> visited_;
    std::vector<uint32_t> unvisited_;
    std::vector<uint32_t> index_in_unvisited_;
    std::vector<bool> is_visited_;
};

#endif /* ifndef ANT_H */
