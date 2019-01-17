import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, ndimage, interpolate, optimize

from mantid import mtd
from mantid.simpleapi import \
    CalculateEfficiencyCorrection, CloneWorkspace, ConvertUnits, \
    ConvertToDistribution, ConvertToHistogram, ConvertToPointData, \
    CreateWorkspace, Divide, LoadNexusMonitors, \
    LoadAscii, Multiply, Rebin, SaveAscii, SplineSmoothing

fitTypeOpts = ['CubicSpline',
               'HowellsFunction',
               'GaussConvCubicSpline',
               'CubicSplineViaMantid']


# Functions for fitting the incident spectrum


def getFitRange(x, y, x_lo, x_hi):
    if x_lo is None:
        x_lo = min(x)
    if x_hi is None:
        x_hi = max(x)

    x_fit = x[(x >= x_lo) & (x <= x_hi)]
    y_fit = y[(x >= x_lo) & (x <= x_hi)]
    return x_fit, y_fit


def fitCubicSpline(x_fit, y_fit, x, s=1e15):
    tck = interpolate.splrep(x_fit, y_fit, s=s)
    fit = interpolate.splev(x, tck, der=0)
    fit_prime = interpolate.splev(x, tck, der=1)
    return fit, fit_prime


def fitCubicSplineViaMantidSplineSmoothing(InputWorkspace,
                                           ParamsInput, ParamsOutput, **kwargs):
    Rebin(
        InputWorkspace=InputWorkspace,
        OutputWorkspace='fit',
        Params=ParamsInput,
        PreserveEvents=True)
    SplineSmoothing(
        InputWorkspace='fit',
        OutputWorkspace='fit',
        OutputWorkspaceDeriv='fit_prime',
        DerivOrder=1,
        **kwargs)
    Rebin(
        InputWorkspace='fit',
        OutputWorkspace='fit',
        Params=ParamsOutput,
        PreserveEvents=True)
    Rebin(
        InputWorkspace='fit_prime_1',
        OutputWorkspace='fit_prime_1',
        Params=ParamsOutput,
        PreserveEvents=True)
    return mtd['fit'].readY(0), mtd['fit_prime_1'].readY(0)


def fitHowellsFunction(x_fit, y_fit, x):
    # Fit with analytical function from HowellsEtAl
    def calc_HowellsFunction(
            lambdas,
            phi_max,
            phi_epi,
            lam_t,
            lam_1,
            lam_2,
            a):
        term1 = phi_max * ((lam_t**4.) / lambdas**5.) * \
            np.exp(-(lam_t / lambdas)**2.)
        term2 = (phi_epi / (lambdas**(1. + 2. * a))) * \
            (1. / (1 + np.exp((lambdas - lam_1) / lam_2)))
        return term1 + term2

    def calc_HowellsFunction1stDerivative(
            lambdas, phi_max, phi_epi, lam_t, lam_1, lam_2, a):
        term1 = (((2 * lam_t**2) / lambdas**2) - 5.) * (1. / lambdas) * \
            phi_max * ((lam_t**4.) / lambdas**5.) * np.exp(-(lam_t / lambdas)**2.)
        term2 = ((1 + 2 * a) / lambdas) \
            * (1. / lambdas) * (phi_epi / (lambdas ** (1. + 2. * a))) \
            * (1. / (1 + np.exp((lambdas - lam_1) / lam_2)))
        return term1 + term2

    params = [1., 1., 1., 0., 1., 1.]
    params, convergence = optimize.curve_fit(
        calc_HowellsFunction, x_fit, y_fit, params)
    fit = calc_HowellsFunction(x, *params)
    fit_prime = calc_HowellsFunction1stDerivative(x, *params)
    return fit, fit_prime


def fitCubicSplineWithGaussConv(x_fit, y_fit, x, sigma=3):
    # Fit with Cubic Spline using a Gaussian Convolution to get weights
    def moving_average(y, sigma=sigma):
        b = signal.gaussian(39, sigma)
        average = ndimage.filters.convolve1d(y, b / b.sum())
        var = ndimage.filters.convolve1d(np.power(y - average, 2), b / b.sum())
        return average, var

    avg, var = moving_average(y_fit)
    spline_fit = interpolate.UnivariateSpline(x_fit, y_fit, w=1. / np.sqrt(var))
    spline_fit_prime = spline_fit.derivative()
    fit = spline_fit(x)
    fit_prime = spline_fit_prime(x)
    return fit, fit_prime


# Get incident spectrum from Monitor

def plotIncidentSpectrum(ax, x, y, x_fit, fit, fit_prime, title=None):
    ax.plot(x, y, 'bo', x_fit, fit, '--', label=title)
    '''
    plt.legend(['Incident Spectrum', 'Fit f(x)'], loc='best')
    if title is not None:
        ax.title(title)
    plt.show()

    plt.plot(x_fit, fit_prime / fit, 'x--', label="Fit f'(x)/f(x)")
    plt.xlabel('Wavelength')
    plt.legend()
    if title is not None:
        plt.title(title)
    axes = plt.gca()
    axes.set_ylim([-12, 6])
    plt.show()
    '''


def FitIncidentSpectrum(InputWorkspace, OutputWorkspace,
                        FitSpectrumWith='GaussConvCubicSpline',
                        BinningForFit="0.15,0.05,3.2",
                        BinningForCalc=None,
                        PlotDiagnostics=False,
                        PlotIncident=False,
                        axis=None):

    if PlotDiagnostics and not axis:
        raise Exception("Must specify a plot axis for plotting diagnostics")

    incident_ws = mtd[InputWorkspace]

    # Fit Incident Spectrum
    # Get axis for actual calc (either provided in BinningForCalc or extracted
    # from incident wksp)
    incident_index = 0
    if BinningForCalc is None:
        x = incident_ws.readX(incident_index)
        binning = [str(i) for i in [min(x), x[1] - x[0], max(x) + x[1] - x[0]]]
        BinningForCalc = ",".join(binning)
    else:
        try:
            params = [float(x) for x in BinningForCalc.split(',')]
        except AttributeError:
            params = [float(x) for x in BinningForCalc]
        xlo, binsize, xhi = params
        x = np.arange(xlo, xhi, binsize)

    Rebin(
        incident_ws,
        OutputWorkspace='fit',
        Params=BinningForFit,
        PreserveEvents=True)
    x_fit = np.array(mtd['fit'].readX(incident_index))
    y_fit = np.array(mtd['fit'].readY(incident_index))

    if len(x_fit) != len(y_fit):
        x_fit = x_fit[:-1]

    if PlotDiagnostics and PlotIncident:
        axis.plot(x_fit, y_fit, 'bo', label='Incident')

    if FitSpectrumWith == 'CubicSpline':
        fit, fit_prime = fitCubicSpline(x_fit, y_fit, x, s=1e7)
        if PlotDiagnostics:
            axis.plot(x, fit, '--', label='Simple Cubic Spline: Default')
            # plotIncidentSpectrum(axis, x_fit, y_fit, x, fit, fit_prime,
            #    title='Simple Cubic Spline: Default')

    elif FitSpectrumWith == 'CubicSplineViaMantid':
        fit, fit_prime = fitCubicSplineViaMantidSplineSmoothing(
            InputWorkspace, ParamsInput=BinningForFit,
            ParamsOutput=BinningForCalc, Error=0.0001, MaxNumberOfBreaks=0)
        if PlotDiagnostics:
            axis.plot(x, fit, '--', label='Cubic Spline via Mantid SplineSmoothing')
            # plotIncidentSpectrum(axis, x_fit, y_fit, x, fit, fit_prime,
            #    title='Cubic Spline via Mantid SplineSmoothing')

    elif FitSpectrumWith == 'HowellsFunction':
        fit, fit_prime = fitHowellsFunction(x_fit, y_fit, x)
        if PlotDiagnostics:
            axis.plot(x, fit, '--', label='HowellsFunction')
            # plotIncidentSpectrum(axis, x_fit, y_fit, x, fit, fit_prime,
            #    title='HowellsFunction')

    elif FitSpectrumWith == 'GaussConvCubicSpline':
        fit, fit_prime = fitCubicSplineWithGaussConv(x_fit, y_fit, x, sigma=0.5)
        if PlotDiagnostics:
            axis.plot(x, fit, '--', label='Cubic Spline w/ Gaussian Kernel Convolution ')
            # plotIncidentSpectrum(axis, x_fit, y_fit, x, fit, fit_prime,
            #    title='Cubic Spline w/ Gaussian Kernel Convolution ')

    else:
        raise Exception("Unknown method for fitting incident spectrum")
        return

    CreateWorkspace(
        DataX=x,
        DataY=np.append(
            fit,
            fit_prime),
        OutputWorkspace=OutputWorkspace,
        UnitX='Wavelength',
        NSpec=2,
        Distribution=False)
    return mtd[OutputWorkspace]


def getIncidentSpectrumParameters():
    incident_spectrums = dict()
    incident_spectrums['Ambient 300K polyethylene'] = {
        'phi_max': 6324,
        'phi_epi': 786,
        'lam_t': 1.58,
        'alpha': 0.099,
        'lam_1': 0.67143,
        'lam_2': 0.06075,
        'filename': 'milder_moderator_polyethlyene_300K'
    }

    incident_spectrums['Ambient 300K poisoned (Gd "meat" in polyethylene slabs)'] = {
        'phi_max': 1200,
        'phi_epi': 786,
        'lam_t': 1.58,
        'alpha': 0.099,
        'lam_1': 0.67143,
        'lam_2': 0.06075,
        'filename': 'milder_moderator_polyethlyene_300K_slabs_gd'
    }

    incident_spectrums['Cold 77K polyethylene (Howells)'] = {
        'phi_max': 3838,
        'phi_epi': 1029,
        'lam_t': 2.97,
        'alpha': 0.089,
        'lam_1': 1.3287,
        'lam_2': 0.14735,
        'filename': 'howells_moderator_polyethlyene_77K'
    }

    incident_spectrums['Cold 77K polyethylene (MildnerEtAl)'] = {
        'phi_max': 3838,
        'phi_epi': 1029,
        'lam_t': 2.97,
        'alpha': 0.089,
        'lam_1': 1.3287,
        'lam_2': 0.14735,
        'a': 0.2882,
        'lam_3': 0.7253,
        'lam_4': 0.0486,
        'filename': 'mildner_moderator_cold_polyethlyene_77K'
    }
    return incident_spectrums


def runIncidentSpectrumTest(plot=False):
    # Howells function
    def delta(lam, lam_1, lam_2, lam_3=None, lam_4=None, alpha=None):
        retVal = 1. / (1. + np.exp((lam - lam_1) / lam_2))

        if lam_3 and lam_4 and alpha:
            term_2 = 1. + (alpha / (1. + np.exp((lam_3 - lam) / lam_4)))
            retVal *= term_2

        return retVal

    def phi_m(lam, **kwargs):
        phi_max = kwargs['phi_max']
        lam_t = kwargs['lam_t']
        return phi_max * (lam_t**4. / lam**5.) * np.exp(-(lam_t / lam)**2.)

    def phi_e(lam, **kwargs):
        phi_epi = kwargs['phi_epi']
        alpha = kwargs['alpha']
        lam_1 = kwargs['lam_1']
        lam_2 = kwargs['lam_2']

        a = None
        lam_3 = None
        lam_4 = None
        if 'a' in kwargs:
            a = kwargs['a']
        if 'lam_3' in kwargs:
            lam_3 = kwargs['lam_3']
        if 'lam_4' in kwargs:
            lam_4 = kwargs['lam_4']

        if lam_3 and lam_4:
            delta_term = delta(lam, lam_1, lam_2, lam_3, lam_4, a)
        else:
            delta_term = delta(lam, lam_1, lam_2)
        return phi_epi * delta_term / (lam**(1 + 2 * alpha))

    def calc_HowellsFunction(lam, **kwargs):
        return phi_m(lam, **kwargs) + phi_e(lam, **kwargs)

    incident_spectrums = getIncidentSpectrumParameters()

    if plot:
        fig, (ax_bm, ax_eff, ax_reload) = plt.subplots(
            3, subplot_kw={'projection': 'mantid'}, sharex=True)

    # Using Mantid
    lam_lo = 0.1
    lam_hi = 3.0
    lam_delta = 0.01
    binning = "%s,%s,%s" % (lam_lo, lam_delta, lam_hi)
    for moderator, spectrum_params in incident_spectrums.items():
        if plot:
            color = next(ax_bm._get_lines.prop_cycler)['color']

        # Create the already efficiency corrected incident spectrum
        incident_ws = 'howells_%s' % moderator
        CreateWorkspace(OutputWorkspace=incident_ws, NSpec=1, DataX=[0], DataY=[0],
                        UnitX='Wavelength', VerticalAxisUnit='Text',
                        VerticalAxisValues='IncidentSpectrum',
                        Distribution=False)
        Rebin(InputWorkspace=incident_ws, OutputWorkspace=incident_ws, Params=binning)
        ConvertToPointData(InputWorkspace=incident_ws, OutputWorkspace=incident_ws)

        wavelengths = mtd[incident_ws].readX(0)
        incident_spectrum = calc_HowellsFunction(wavelengths, **spectrum_params)
        mtd[incident_ws].setY(0, incident_spectrum)

        if plot:
            ax_bm.plot(mtd[incident_ws], '-', color=color, wkspIndex=0, label=moderator)

        # Calculate the efficiency correction and back-calculate the measurement
        # before efficiency correction
        eff_ws = 'efficiency'
        CalculateEfficiencyCorrection(InputWorkspace=incident_ws,
                                      Alpha=0.693,
                                      OutputWorkspace=eff_ws)
        ConvertToPointData(InputWorkspace=eff_ws, OutputWorkspace=eff_ws)

        if plot:
            ax_eff.plot(mtd[eff_ws], '-', color=color, wkspIndex=0, label=moderator + ' efficiency')

        sample_ws = 'sample_ws'
        Divide(LHSWorkspace=incident_ws, RHSWorkspace=eff_ws, OutputWorkspace=sample_ws)

        if plot:
            ax_bm.plot(mtd[sample_ws], 'o', color=color, wkspIndex=0,
                       label=moderator + ' measurement')

        # Save the initial, pre-corrected spectrum, load back in and then re-apply
        # the correction for test
        reloaded_ws = 'reloaded_ws'
        reloaded_corr_ws = 'reloaded_corr_ws'
        savefile = spectrum_params['filename'] + '.txt'
        SaveAscii(InputWorkspace=sample_ws, Filename=savefile)

        LoadAscii(Filename=savefile, OutputWorkspace=reloaded_ws, Unit='Wavelength')

        CalculateEfficiencyCorrection(InputWorkspace=reloaded_ws,
                                      Alpha=0.693,
                                      OutputWorkspace=eff_ws)
        Multiply(LHSWorkspace=reloaded_ws,
                 RHSWorkspace=eff_ws,
                 OutputWorkspace=reloaded_corr_ws)

        if plot:
            ax_reload.plot(mtd[reloaded_ws], 'o', color=color, wkspIndex=0,
                           label=moderator + 'reloaded measurement')
            ax_reload.plot(mtd[reloaded_corr_ws], '-', color=color, wkspIndex=0,
                           label=moderator + 'reloaded corrected')

    if plot:
        ax_eff.legend()
        for axis in [ax_bm, ax_reload]:
            axis.legend()
            axis.set_ylim(0.0, 6000.0)
        plt.show()


def runNomadTest(plot=False):
    config = dict()
    config['Instrument'] = "NOM"
    config['Sample'] = {"Runs": "33943",
                        "InelasticCorrection": {"Type": "Placzek",
                                                "Order": "1st",
                                                "Self": True,
                                                "Interference": False,
                                                "FitSpectrumWith": "GaussConvCubicSpline",
                                                "LambdaBinningForFit": "0.16,0.04,2.8",
                                                "LambdaBinningForCalc": "0.1,0.0001,3.0"}
                        }
    binning = "0.0212406,0.00022,3.39828"  # matches read_bm.pro for lambda[100:15999]
    sample = config['Sample']

    # Get incident spectrum test for NOMAD
    runs = sample["Runs"].split(',')
    runs = ["%s_%s" % (config["Instrument"], run) for run in runs]

    monitor = 'monitor'
    incident_ws = 'incident_ws'

    if plot:
        fig, (ax_bm, ax_bmeff) = plt.subplots(2, subplot_kw={'projection': 'mantid'}, sharex=True)

    # Beam Monitor
    LoadNexusMonitors(Filename=runs[0], OutputWorkspace=monitor)
    ConvertUnits(InputWorkspace=monitor, OutputWorkspace=monitor,
                 Target='Wavelength', EMode='Elastic')

    Rebin(InputWorkspace=monitor,
          OutputWorkspace=monitor,
          Params=binning,
          PreserveEvents=False)

    ConvertToPointData(InputWorkspace=monitor, OutputWorkspace=monitor)

    if plot:
        ax_bm.plot(mtd[monitor], '-', wkspIndex=0, label='Monitor')
        ax_bm.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    # Use sample info
    CalculateEfficiencyCorrection(WavelengthRange=binning,
                                  ChemicalFormula="(He3)",
                                  DensityType="Number Density",
                                  Density=1.93138101e-08,
                                  Thickness=.1,
                                  OutputWorkspace=incident_ws)

    Multiply(LHSWorkspace=monitor, RHSWorkspace=incident_ws, OutputWorkspace=incident_ws)
    mtd[incident_ws].setDistribution(True)

    if plot:
        ax_bmeff.plot(mtd[incident_ws], '-', wkspIndex=0, label='Incident Spectrum (density)')

    # Use measured efficiency
    CalculateEfficiencyCorrection(WavelengthRange=binning,
                                  ChemicalFormula="(He3)",
                                  MeasuredEfficiency=1.03e-5,
                                  OutputWorkspace=incident_ws)

    Multiply(LHSWorkspace=monitor, RHSWorkspace=incident_ws, OutputWorkspace=incident_ws)
    mtd[incident_ws].setDistribution(True)

    if plot:
        ax_bmeff.plot(mtd[incident_ws], 'o', wkspIndex=0, label='Incident Spectrum (efficiency)')

    # Use alpha
    CalculateEfficiencyCorrection(WavelengthRange=binning,
                                  Alpha=5.72861786781e-06,
                                  OutputWorkspace=incident_ws)

    Multiply(LHSWorkspace=monitor, RHSWorkspace=incident_ws, OutputWorkspace=incident_ws)
    mtd[incident_ws].setDistribution(True)

    if plot:
        ax_bmeff.plot(mtd[incident_ws], '--', wkspIndex=0, label='Incident Spectrum (alpha)')

    '''
    # Rebin workspace
    Rebin(InputWorkspace=incident_ws,
          OutputWorkspace=incident_ws,
          Params="0.2,0.01,4.0",
          PreserveEvents=False)

    if plot:
        ax_bmeff.plot(mtd[incident_ws], '--', wkspIndex=0, label='Incident Spectrum (alpha + rebin')
    '''

    # Plot all
    if plot:
        ax_bm.legend()
        ax_bmeff.legend()
        plt.show()


def runFitIncidentSpectrumTest(axis=None, axis_indices=None, plot=False, show_plot=False):
    # Fit incident spectrum
    incident_spectrums = getIncidentSpectrumParameters()
    incident_ws = 'incident_ws'
    eff_ws = 'efficiency'
    savefile = incident_spectrums['Ambient 300K polyethylene']['filename']
    LoadAscii(Filename=savefile, OutputWorkspace=incident_ws, Unit='Wavelength')
    CalculateEfficiencyCorrection(InputWorkspace=incident_ws,
                                  Alpha=0.693,
                                  OutputWorkspace=eff_ws)
    ConvertToPointData(InputWorkspace=incident_ws,
                       OutputWorkspace=incident_ws)
    Multiply(LHSWorkspace=incident_ws,
             RHSWorkspace=eff_ws,
             OutputWorkspace=incident_ws)
    mtd[incident_ws].setDistribution(True)

    incident_fit_prefix = 'incident_fit'
    if plot:
        if axis is None:
            fig, axis = plt.subplots(5, subplot_kw={'projection': 'mantid'},
                                     sharex=True, sharey=True)
    if axis_indices is None:
        if axis is None:
            axis_indices = range(5)
        else:
            axis_indices = range(len(axis))

    rebin = True
    if plot:
        first_axis = axis[axis_indices[0]]
    else:
        first_axis = None

    plotIncident = True
    for idx, fit_type in zip(axis_indices[1:], fitTypeOpts):

        incident_fit = incident_fit_prefix + "_" + fit_type
        FitIncidentSpectrum(InputWorkspace=incident_ws,
                            OutputWorkspace=incident_fit,
                            FitSpectrumWith=fit_type,
                            BinningForFit="0.1,0.01,3.0",
                            BinningForCalc="0.1,0.0001,3.0",
                            PlotDiagnostics=plot,
                            PlotIncident=plotIncident,
                            axis=first_axis)
        plotIncident = False

        incident_ws_rebin = "incident_ws_rebin"
        if rebin:
            Rebin(incident_ws,
                  OutputWorkspace=incident_ws_rebin,
                  Params="0.15,0.05,3.2",
                  PreserveEvents=True)
        else:
            CloneWorkspace(InputWorkspace=incident_ws,
                           OutputWorkspace=incident_ws_rebin)
            ConvertToHistogram(InputWorkspace=incident_ws_rebin,
                               OutputWorkspace=incident_ws_rebin)
            ConvertToDistribution(Workspace=incident_ws_rebin)

        if plot:
            axis[idx].plot(mtd[incident_ws], '-', wkspIndex=0, label="Incident")
            axis[idx].plot(mtd[incident_ws_rebin], '-', wkspIndex=0, label="Incident Rebinned")
            axis[idx].plot(mtd[incident_fit], '-', wkspIndex=0, label=fit_type)
            axis[idx].legend()
            axis[axis_indices[0]].legend()

    if show_plot:
        plt.show()


def runFitNomadIncidentSpectrumTest(axis=None, axis_indices=None, plot=False, show_plot=False):
    incident_ws = 'incident_ws'
    incident_fit = 'incident_fit'
    if plot:
        if axis is None:
            fig, axis = plt.subplots(
                5, subplot_kw={
                    'projection': 'mantid'}, sharex=True, sharey=True)
    if axis_indices is None:
        if axis is None:
            axis_indices = range(5)
        else:
            axis_indices = range(len(axis))

    rebin = True
    if plot:
        first_axis = axis[axis_indices[0]]
    else:
        first_axis = None

    plotIncident = True
    for idx, fit_type in zip(axis_indices[1:], fitTypeOpts):

        FitIncidentSpectrum(InputWorkspace=incident_ws,
                            OutputWorkspace=incident_fit,
                            FitSpectrumWith=fit_type,
                            BinningForFit="0.02,0.01,3.0",
                            PlotDiagnostics=plot,
                            PlotIncident=plotIncident,
                            axis=first_axis)
        plotIncident = False
        # BinningForFit="0.16,0.04,2.8",
        # BinningForCalc="0.1,0.0001,3.0",

        incident_ws_rebin = "incident_ws_rebin"
        if rebin:
            Rebin(incident_ws,
                  OutputWorkspace=incident_ws_rebin,
                  Params="0.0212406,0.00022,3.39828",  # matches read_bm.pro for lambda[100:15999]
                  PreserveEvents=True)
        else:
            CloneWorkspace(InputWorkspace=incident_ws,
                           OutputWorkspace=incident_ws_rebin)
            ConvertToHistogram(InputWorkspace=incident_ws_rebin,
                               OutputWorkspace=incident_ws_rebin)
            ConvertToDistribution(Workspace=incident_ws_rebin)

        if plot:
            axis[idx].plot(mtd[incident_ws], '-', wkspIndex=0, label="Incident")
            axis[idx].plot(mtd[incident_ws_rebin], 'x', wkspIndex=0, label="Incident Rebinned")
            axis[idx].plot(mtd[incident_fit], '-', wkspIndex=0, label=fit_type)
            axis[idx].legend()
            axis[axis_indices[0]].legend()
            axis[idx].set_ylim(0.0, 2e9)

    if show_plot:
        plt.show()


if '__main__' == __name__:
    howellsTestFlag = True
    nomadTestFlag = False
    fitIncidentSpectrumFlag = True
    fitNomadIncidentSpectrumFlag = False
    plotFlag = True

    # Howells incident spectrum for testing
    if howellsTestFlag:
        runIncidentSpectrumTest(plot=plotFlag)

    # Nomad incident spectrum from monitor
    if nomadTestFlag:
        runNomadTest(plot=plotFlag)

    # Fitting the Howells spectrums
    if fitIncidentSpectrumFlag:
        fig, ax = plt.subplots(3, 2, subplot_kw={'projection': 'mantid'}, sharex=True, sharey=True)
        idx = [(0, 0), (1, 0), (1, 1), (2, 0), (2, 1)]
        runIncidentSpectrumTest()
        runFitIncidentSpectrumTest(axis=ax, axis_indices=idx, plot=plotFlag, show_plot=plotFlag)

    # Fitting the Nomad monitor
    if fitNomadIncidentSpectrumFlag:
        fig, ax = plt.subplots(3, 2, subplot_kw={'projection': 'mantid'}, sharex=True, sharey=True)
        idx = [(0, 0), (1, 0), (1, 1), (2, 0), (2, 1)]
        runNomadTest()
        runFitNomadIncidentSpectrumTest(axis=ax, axis_indices=idx,
                                        plot=plotFlag, show_plot=plotFlag)
