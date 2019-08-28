
def set_inelastic_correction(inelastic_dict):
    default_inelastic_dict = {"Type": None}

    if inelastic_dict is None:
        return default_inelastic_dict

    corr_type = inelastic_dict["Type"]
    if corr_type is None or corr_type == u'None':
        return default_inelastic_dict

    if corr_type:
        if corr_type == "Placzek":
            default_settings = {"Order": "1st",
                                "Self": True,
                                "Interference": False,
                                "FitSpectrumWith": "GaussConvCubicSpline",
                                "LambdaBinning": "0.16,0.04,2.8"}
            inelastic_settings = default_settings.copy()
            inelastic_settings.update(inelastic_dict)

        else:
            raise Exception("Unknown Inelastic Correction Type")

    return inelastic_settings
    