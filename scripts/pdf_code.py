    SNSPowderReduction(
        Filename=sam_scans,
        MaxChunkSize=alignAndFocusArgs['MaxChunkSize'],
        PreserveEvents=True,
        PushDataPositive='ResetToZero',
        CalibrationFile=alignAndFocusArgs['CalFilename'],
        CharacterizationRunsFile=merging['Characterizations']['Filename'],
        BackgroundNumber=sample["Background"]["Runs"],
        VanadiumNumber=van["Runs"],
        VanadiumBackgroundNumber=van["Background"]["Runs"],
        RemovePromptPulseWidth=alignAndFocusArgs['RemovePromptPulseWidth'],
        ResampleX=alignAndFocusArgs['ResampleX'],
        BinInDspace=True,
        FilterBadPulses=25.,
        SaveAs="gsas fullprof topas",
        OutputFilePrefix=title,
        OutputDirectory=OutputDir,
        StripVanadiumPeaks=True,
        VanadiumRadius=van_geometry['Radius'],
        NormalizeByCurrent=True,
        FinalDataUnits="dSpacing")

    # Ouput bank-by-bank with linear fits for high-Q

    # fit the last 80% of the bank being used
    for i, q in zip(range(mtd[sam_corrected].getNumberHistograms()), qmax):
        qmax_data = getQmaxFromData(sam_corrected, i)
        qmax[i] = q if q <= qmax_data else qmax_data

    fitrange_individual = [(high_q_linear_fit_range*q, q) for q in qmax]

    for q in qmax:
        print('Linear Fit Qrange:', high_q_linear_fit_range*q, q)


    kwargs = { 'btot_sqrd_avg' : btot_sqrd_avg,
               'bcoh_avg_sqrd' : bcoh_avg_sqrd,
               'self_scat' : self_scat }

    save_banks_with_fit( title, fitrange_individual, InputWorkspace='SQ_banks', **kwargs)
    save_banks_with_fit( title, fitrange_individual, InputWorkspace='FQ_banks', **kwargs)
    save_banks_with_fit( title, fitrange_individual, InputWorkspace='FQ_banks_raw', **kwargs)

    save_banks('SQ_banks', title=os.path.join(OutputDir,title+"_SQ_banks.dat"),
                binning=binning)
    save_banks('FQ_banks', title=os.path.join(OutputDir,title+"_FQ_banks.dat"),
                binning=binning)
    save_banks('FQ_banks_raw', title=os.path.join(OutputDir,title+"_FQ_banks_raw.dat"),
                binning=binning)

    # Event workspace -> Histograms
    Rebin(InputWorkspace=sam_corrected, OutputWorkspace=sam_corrected,
          Params=binning, PreserveEvents=True)
    Rebin(InputWorkspace=van_corrected, OutputWorkspace=van_corrected,
          Params=binning, PreserveEvents=True)
    Rebin(InputWorkspace='container',   OutputWorkspace='container',
          Params=binning, PreserveEvents=True)
    Rebin(InputWorkspace='sample', OutputWorkspace='sample',
          Params=binning, PreserveEvents=True)
    if van_bg is not None:
        Rebin(InputWorkspace=van_bg, OutputWorkspace='background',
              Params=binning, PreserveEvents=True)

    # Apply Qmin Qmax limits

    #MaskBinsFromTable(InputWorkspace=sam_corrected, OutputWorkspace='sam_single',
                       MaskingInformation=mask_info)
    #MaskBinsFromTable(InputWorkspace=van_corrected, OutputWorkspace='van_single',
                       MaskingInformation=mask_info)
    #MaskBinsFromTable(InputWorkspace='container', OutputWorkspace='container_single',
                       MaskingInformation=mask_info)
    #MaskBinsFromTable(InputWorkspace='sample', OutputWorkspace='sample_raw_single',
                       MaskingInformation=mask_info)

    # Get sinlge, merged spectrum from banks

    CloneWorkspace(InputWorkspace=sam_corrected, OutputWorkspace='sam_single')
    CloneWorkspace(InputWorkspace=van_corrected, OutputWorkspace='van_single')
    CloneWorkspace(InputWorkspace='container', OutputWorkspace='container_single')
    CloneWorkspace(InputWorkspace='sample', OutputWorkspace='sample_raw_single')
    CloneWorkspace(InputWorkspace='background', OutputWorkspace='background_single')

    SumSpectra(InputWorkspace='sam_single', OutputWorkspace='sam_single',
               ListOfWorkspaceIndices=wkspIndices)
    SumSpectra(InputWorkspace='van_single', OutputWorkspace='van_single',
               ListOfWorkspaceIndices=wkspIndices)

    # Diagnostic workspaces
    SumSpectra(InputWorkspace='container_single', OutputWorkspace='container_single',
               ListOfWorkspaceIndices=wkspIndices)
    SumSpectra(InputWorkspace='sample_raw_single', OutputWorkspace='sample_raw_single',
               ListOfWorkspaceIndices=wkspIndices)
    SumSpectra(InputWorkspace='background_single', OutputWorkspace='background_single',
               ListOfWorkspaceIndices=wkspIndices)

    # Merged S(Q) and F(Q)

    save_banks(InputWorkspace="FQ_banks_ws",
               Filename=nexus_filename,
               Title="FQ_banks",
               OutputDir=OutputDir,
               Binning=binning)
    save_banks(InputWorkspace="SQ_banks_ws",
               Filename=nexus_filename,
               Title="SQ_banks",
               OutputDir=OutputDir,
               Binning=binning)


    # do the division correctly and subtract off the material specific term
    CloneWorkspace(InputWorkspace='sam_single', OutputWorkspace='SQ_ws')
    SQ = (1./bcoh_avg_sqrd)*mtd['SQ_ws'] - (term_to_subtract-1.)  # +1 to get back to S(Q)

    CloneWorkspace(InputWorkspace='sam_single', OutputWorkspace='FQ_ws')
    FQ_raw = mtd['FQ_ws']
    FQ = FQ_raw - self_scat

    qmax = 48.0
    Fit(Function='name=LinearBackground,A0=1.0,A1=0.0',
        StartX=high_q_linear_fit_range*qmax, EndX=qmax, # range cannot include area with NAN
        InputWorkspace='SQ', Output='SQ', OutputCompositeMembers=True)
    fitParams = mtd['SQ_Parameters']

    qmax = getQmaxFromData('FQ', WorkspaceIndex=0)
    Fit(Function='name=LinearBackground,A0=1.0,A1=0.0',
        StartX=high_q_linear_fit_range*qmax, EndX=qmax, # range cannot include area with NAN
        InputWorkspace='FQ', Output='FQ', OutputCompositeMembers=True)
    fitParams = mtd['FQ_Parameters']

    qmax = 48.0
    Fit(Function='name=LinearBackground,A0=1.0,A1=0.0',
        StartX=high_q_linear_fit_range*qmax, EndX=qmax, # range cannot include area with NAN
        InputWorkspace='FQ_raw', Output='FQ_raw', OutputCompositeMembers=True)
    fitParams = mtd['FQ_raw_Parameters']

    # Save dat file
    header_lines = ['<b^2> : %f ' % btot_sqrd_avg, \
                    '<b>^2 : %f ' % bcoh_avg_sqrd, \
                    'self scattering: %f ' % self_scat, \
                    'fitrange: %f %f '  % (high_q_linear_fit_range*qmax,qmax), \
                    'for merged banks %s: %f + %f * Q' \
                    % (','.join([ str(i) for i in wkspIndices]), \
                    fitParams.cell('Value', 0), fitParams.cell('Value', 1)) ]
