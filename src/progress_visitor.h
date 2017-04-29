#ifndef PROGRESS_VISITOR_H
#define PROGRESS_VISITOR_H

/*
 * A base class for all algorithm progress trackers
 */
template<typename Alg>
class ProgressVisitor {
public:

    virtual ~ProgressVisitor() {}

    /*
     * This method is called after the algorithm was initialized and ready to
     * perform first iteration.
     */
    virtual void started(Alg &alg) = 0;


    /*
     * This method is called after the algorithm has finished next iteration.
     */
    virtual void iteration_done(Alg &alg) = 0;


    /*
     * This method is called after all the iterations are completed.
     */
    virtual void finished(Alg &alg) = 0;
};

#endif /* ifndef PROGRESS_VISITOR_h */
