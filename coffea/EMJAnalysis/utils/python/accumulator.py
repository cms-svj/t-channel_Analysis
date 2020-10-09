from coffea import hist
import coffea.processor as processor

def print_dict_accumulator(d, keys=None):
    print(type(d))
    print(d)
    if keys != None:
        for key in keys:
            print(d[key])
            print(d[key].values())
            print(d[key].axis('dataset').identifiers())
            print(d[key].identifiers('dataset'))
    print('\n')   

def filter_histograms(histograms, histogram_names=None, search=None, histograms_only=False):
    output = histograms
    if histograms_only:
        output = processor.accumulator.dict_accumulator(filter(lambda x: type(x[1]) == hist.hist_tools.Hist, output.items()))    
    if histogram_names == None and search == None:
        return output
    else:
        if histogram_names != None:
            output = processor.accumulator.dict_accumulator({ hname: output[hname] for hname in histogram_names if hname in output })
        if search != None:
            # Get the first histogram in the dictionary
            # This could break if the first entry isn't a histogram
            first_hist_key = list(output.keys())[0]
            identifiers = output[first_hist_key].identifiers('dataset')
            grouping = { identifier.name: [identifier.name] for identifier in identifiers if any(search_string in identifier.name for search_string in search) }
            output = processor.accumulator.dict_accumulator({ key:(value.group('dataset', hist.Cat('dataset', 'Primary dataset', sorting='placement'), grouping) if type(value) == hist.hist_tools.Hist else value) for (key, value) in output.items() })
            # You can also do:
            #    output_new_single_hist = output['<hist_name>'][idnew,:]
            # Where idnew is a list of process named (identifiers)
            # Note, this will not filter the identifiers on the 'dataset' axis, just in the available values
    return output

def test_filter_histograms(accumulator_dict):
    selected = filter_histograms(accumulator_dict)
    print_dict_accumulator(selected,['jtpt'])

    selected = filter_histograms(accumulator_dict,histograms_only=True)
    print_dict_accumulator(selected,['jtpt'])

    selected = filter_histograms(accumulator_dict,histogram_names=['jtpt'])
    print_dict_accumulator(selected,['jtpt'])

    selected = filter_histograms(accumulator_dict,search=['ctau'])
    print_dict_accumulator(selected,['jtpt','jteta'])

    selected = filter_histograms(accumulator_dict,histogram_names=['jtpt','jteta','cutflow'],search=['ctau'],histograms_only=True)
    print_dict_accumulator(selected,['jtpt','jteta'])
    
def integrate_hist_over_dataset(accumulator_dict, histogram_name, processes):
    return accumulator_dict[histogram_name].integrate('dataset',processes)

def test_integrate_hist_over_dataset(accumulator_dict):
    tmp = integrate_hist_over_dataset(accumulator_dict,'jtpt',["EMJ_2016_mMed-1000_mDark-20_ctau-1000_unflavored-down"])
    print(tmp.values())
    tmp = integrate_hist_over_dataset(accumulator_dict,'jtpt',["EMJ_2016_mMed-1000_mDark-20_kappa-0p12_aligned-down"])
    print(tmp.values())
    tmp = integrate_hist_over_dataset(accumulator_dict,'jtpt',["EMJ_2016_mMed-1000_mDark-20_ctau-1000_unflavored-down","EMJ_2016_mMed-1000_mDark-20_kappa-0p12_aligned-down"])
    print(tmp.values())