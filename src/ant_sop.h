#ifndef ANT_SOP_H
#define ANT_SOP_H

#include "sop_instance.h"
#include "ant.h"

/*
 * Ant class adjusted to the SOP problem.
 *
 * While building a solution to the SOP we need to consider precedence
 * constraints which are not present in the TSP.
 */
class AntSOP final : public Ant {
public:
    AntSOP(const SOPInstance &problem);

    virtual ~AntSOP() {}

    void move_to(uint32_t node) override;

    bool is_available(uint32_t node) const override { return available_[node]; }

private:

    const SOPInstance &problem_;
    std::vector<uint32_t> incoming_count_; // incoming_count_[v] = number of remaining edges incoming to v, i.e. (u, v)
    std::vector<bool> available_; // available_[v] = 1 if node is not a part of solution (yet), 0 otherwise
};

#endif /* ifndef ANT_SOP_H */
