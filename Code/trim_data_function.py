def data(data, band1, band2, band3, band4):
    '''Select data for analysis'''
    ratio = 0.25
    data = data#[10000:30000]
    # Colour 1
    wave1 = data[band1]
    wave1_unc = data[band1+'_unc']
    wave2 = data[band2]
    wave2_unc = data[band2+'_unc']
    # Colour 2
    wave3 = data[band3]
    wave3_unc = data[band3+'_unc']
    wave4 = data[band4]
    wave4_unc = data[band4+'_unc']
    # Change parameters to match data_file
    # Remove data pieces with no value

    wave1_trim = np.logical_and(wave1 != -99, wave1_unc != -99)
    wave2_trim = np.logical_and(wave2 != -99, wave2_unc != -99)
    wave3_trim = np.logical_and(wave3 != -99, wave3_unc != -99)
    wave4_trim = np.logical_and(wave4 != -99, wave4_unc != -99)

    colour1_ratio = np.logical_and(wave1_unc < ratio,
                                   wave2_unc < ratio)
    colour2_ratio = np.logical_and(wave3_unc < ratio,
                                   wave4_unc < ratio)

    gooddata1 = np.logical_and(np.logical_and(wave1_trim, wave2_trim),
                               np.logical_and(wave3_trim, wave4_trim))

    # Remove data above certain magnitude
    gooddata2 = np.logical_and(colour1_ratio, colour2_ratio)

    # Only data that match criteria for both colours
    greatdata = np.logical_and(gooddata1, gooddata2)

    colour1 = wave1[greatdata] - wave2[greatdata]
    colour2 = wave3[greatdata] - wave4[greatdata]
    
    return(wave1, wave2, wave3, wave4, greatdata, colour1, colour2)

def make_directory(p_path):
    '''Save results of each analysis
            save_path: set path of where you would like results saved'''
    # Create new plots_folder
    pl_path = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\{}'.format(
              p_path)
    if not os.path.exists(pl_path):
        os.makedirs(pl_path)

    return(pl_path)

def load_data(surveyfile_, experiments):
    '''User upload data file'''

    d_file_name = str(surveyfile_)
    t_file_name = str(experiments)
    data = Table.read(d_file_name, format='ascii.commented_header',
                      guess=False)
    tests = Table.read(t_file_name, format='ascii.commented_header',
                       guess=False)

    return (data, tests)