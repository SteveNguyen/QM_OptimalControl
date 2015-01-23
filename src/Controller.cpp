#include "Controller.h"

Controller::Controller(Optimizer& eng_, ModelLearner& learner_) : eng(eng_), learner(learner_), state(eng_.getStateVar()), action(eng_.getActionVar())
{



}

Controller::~Controller()
{
    //dtor
}





