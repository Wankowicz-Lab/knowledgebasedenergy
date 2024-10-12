def get_energy_from_stats(potential, measured, default_maximum=True,
                          set_maximum=None):
    ‘’'
    Arguments:
        potential (pandas dataframe): potential energy dataframe wih columns:
            bin min, bin max and E
        measured (float): measured geometric parameter to get energy of
        default_maximum (bool): if True, maximum eneregy will be assigned to
            geometric parameters not found in any bins that the potential energy
            function covers. If False, will use the value specified in set_maximum {default: True}
        set_maximum (float): set maximum energy value when the geometric
            parameter is not found in any bins that the potential energy function covers.
    Returns
        energy value (kcal/mol).
    ‘’'
    if np.isnan(measured):
        return np.nan
    dfE = potential.loc[(potential[‘bin min’] < measured)
                        & (potential[‘bin max’] >= measured)]
    if len(dfE) == 0:
        if default_maximum:
            return potential[‘E’].max()
        elif (not default_maximum) and set_maximum:
            return set_maximum
        else:
            print(“Need to specify set_maximum if not using default maximum energy”)
    else:
        return dfE.iloc[0][‘E’]
