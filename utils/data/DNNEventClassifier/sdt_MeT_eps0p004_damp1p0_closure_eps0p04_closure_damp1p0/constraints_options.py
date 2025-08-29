#Options for the constraints
#Constraints dictionary
#code assumes x axis is DNN score
#To use closure constraint must also set constraint on y axis, and must be last constraint in list
#here max = maximum inequality constraint, min = minimum inequality constraint, eq = equality constraint

features_to_constrain = ["MET"]
# Number of constraints to use, if do not want to use constraints, then set to 0
constraints = {
    "MET": {
        "axis": "y",
        "type": "max",
        "epsilon": 0.004,
        "damping": 1.0,
        "scale": 1.0,
        "power": 1.0,
        "with_signal": False,
        "reweighted": False,
    },
    "closure": {
        "axis": None,
        "type": "max",
        "epsilon": 0.04,
        "damping": 1.0,
        "scale": 1.0,
        "symmetric": True,
        "reweighted": False,           #here it reweights the actual number of events in the batch to the expected number of events (x-section) [only for closure loss calculation]
        "n_events_min": 1,     # minimum number of MC events in each ABCD region
    },
}
n_contraints = len(constraints.keys())

#BCE is always monitored, e.g. to use when still want to calculate losses which are not constrained (e.g. closure loss, mass regression loss, etc.)
#here don't put losses you want to constrain, put only losses you want to further monitor 
losses_to_monitor = []

