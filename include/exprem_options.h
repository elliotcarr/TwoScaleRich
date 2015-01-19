struct exprem_options_struct
{
    double AbsTol;
    double RelTol;
    int    MinKrylov;
    int    MaxKrylov;
    double MinStep;
    double MaxStep;
    int    MaxNumSteps;
    double SafetyFac;
    double MinStepDecFac;
    double MaxStepIncFac;
    double InitStep;
};