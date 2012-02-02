get_exp_colors <- function() {
    colors = list(control = rgb(0,0,0), control_light = rgb(0,0,0,0.5),
        DMSO = rgb(1,0,0), DMSO_light = rgb(1,0,0,0.5),
        BB94 = rgb(47/255,143/255,51/255), BB94_light = rgb(24/255,82/255,27/255,0.5),
        PP2 = rgb(77/255,0/255,104/255), PP2_light = rgb(77/255,0/255,104/255,0.5),
        Noc = rgb(255/255,128/255,0/255), Noc_light = rgb(255/255,128/255,0/255,0.5),
        PurvA = rgb(0/255,255/255,241/255), PurvA_light = rgb(0/255,255/255,241/255,0.5),
        FAK = rgb(14/255,19/255,149/255), FAK_light = rgb(14/255,19/255,149/255,0.5)
        )

    return(colors)
}
