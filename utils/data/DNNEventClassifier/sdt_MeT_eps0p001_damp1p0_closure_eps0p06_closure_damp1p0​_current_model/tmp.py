import awkward as ak
import numpy as np


def __filter_nans(ak_array):

    filter = None
    for field in ak_array.fields:
        n_dims = str(ak.type(ak_array[field])).count("*")
        if n_dims == 1:
            filter_ = (
                (~np.isnan(ak_array[field]))
                & (~np.isinf(ak_array[field]))
            )
            if filter is None:
                filter = filter_
            else:
                filter = filter & filter_

    new_ak_array = ak_array[filter]

    if len(new_ak_array[field]) != len(ak_array[field]):
        fraction = 1 -  (len(new_ak_array) / len(ak_array))
        print(f"NaN / inf cut removed {fraction}% of the events")

    return new_ak_array
        

ak_array = ak.Array({
    #"a": [[1, 3], [4, 5, 6]],
    #"b": [[1, np.inf], [4, 5, 6]],
    "a": [1, 2, 3, 4, 5, 6],
    "b": [1, 2, np.inf, 4, 5, 6],
})

x = __filter_nans(ak_array)
print(x)

