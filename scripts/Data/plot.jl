using XLSX
using Plots

# Load the Excel file
xlsx_path = "/Users/kneenaugh/Desktop/Git/AstroChemNetwork/plots"  # <- update this path as needed

# Open the Excel file and read the data
XLSX.openxlsx(xlsx_path) do xf

    nel_fid = xf["nelson fiducial"]
    nel_neu = xf["nelson neutral"]

    glo_fid = xf["glover fiducial"]
    glo_neu = xf["glover neutral"]

    umi_fid = xf["umist fiducial"]
    umi_neu = xf["umist neutral"]

    # Read x-axis data (time): A2 to A102
    time = [nel_fid["A$(i)"] for i in 2:335]

    # Read y-axis data (C): G2 to G102
    C_nel_fid = [nel_fid["G$(i)"] for i in 2:335]
    C_nel_neu = [nel_neu["G$(i)"] for i in 2:335]

    C_glo_fid = [glo_fid["K$(i)"] for i in 2:335]
    C_glo_neu = [glo_neu["K$(i)"] for i in 2:335]

    C_umi_fid = [umi_fid["C$(i)"] for i in 2:335]
    C_umi_neu = [umi_neu["C$(i)"] for i in 2:335]

    O_nel_fid = [nel_fid["I$(i)"] for i in 2:335]
    O_nel_neu = [nel_neu["I$(i)"] for i in 2:335]

    O_glo_fid = [glo_fid["M$(i)"] for i in 2:335]
    O_glo_neu = [glo_neu["M$(i)"] for i in 2:335]

    O_umi_fid = [umi_fid["S$(i)"] for i in 2:335]
    O_umi_neu = [umi_neu["S$(i)"] for i in 2:335]

    # color pallettes
    blues = palette(:devon, 10)
    golds = palette(:lajolla10, 10)

    # x-tick labels
    xtick_vals = [7.884e12, 1.5768e13, 2.3652e13, 3.1536e13]
    xtick_labels = ["2.5e5", "5e5", "7.5e5", "1e6"]
    ytick_vals = [0.0, 2.5e-5, 5e-5, 7.5e-5, 1e-4]
    ytick_labels = ["0.0", "2.5e-5", "5e-5", "7.5e-5", "1e-4"]


    ytick_vals2 = [0.0, 5e-5, 1e-4, 1.5e-4]
    ytick_labels2 = ["0.0", "5e-5", "1e-4", "1.5e-4"]




    # Plot the data
    #=
    ### Nelson Vs. Umist: Carbon ###
    plot(time, C_nel_fid, label="Nelson fiducial u0", lw=4, lc = blues[1])
    plot!(time, C_nel_neu, label="Nelson neutral u0", lw=4, lc = blues[4])
    plot!(time, C_umi_fid, label="Umist fiducial u0", lw=4, lc =blues[7])
    plot!(time, C_umi_neu, label="Umist neutral u0", lw=4, lc = blues[9], 
          xlabel="Time (Years)", 
          ylabel="Abundance per Hydrogen Molecule", 
          title="Nelson Vs. Umist: Carbon", 
          legend=:topright, 
          xticks=(xtick_vals, xtick_labels), 
          yticks=(ytick_vals, ytick_labels),
          legendfontsize=10, 
          xguidefontsize=14, 
          yguidefontsize=14, 
          ytickfontsize=10,
          xtickfontsize=10,
          fontfamily="Bookman Light", 
          dpi=300)
    savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/scripts/Data/nelson_umist_carbon_abundance_plot.svg")
=#
    
    #=
    ### Nelson Vs. Umist: Oxygen ###
    plot(time, O_nel_fid, label="Nelson fiducial u0", lw=4, lc = golds[1])
    plot!(time, O_nel_neu, label="Nelson neutral u0", lw=4, lc = golds[4])
    plot!(time, O_umi_fid, label="Umist fiducial u0", lw=4, lc = golds[7])
    plot!(time, O_umi_neu, label="Umist neutral u0", lw=4, lc = golds[9], 
          xlabel="Time (Years)", 
          ylabel="Abundance per Hydrogen Molecule", 
          title="Nelson Vs. Umist: Oxygen", 
          legend=:bottomright, 
          xticks=(xtick_vals, xtick_labels),
          yticks=(ytick_vals2, ytick_labels2),
          legendfontsize=10, 
          xguidefontsize=14, 
          yguidefontsize=14, 
          ytickfontsize=10,
          xtickfontsize=10,
          fontfamily="Bookman Light", 
          dpi=300)
    savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/scripts/Data/nelson_umist_oxygen_abundance_plot.svg")
=#
    #=
    ### Glover Vs. Umist: Carbon ###
    plot(time, C_glo_fid, label="Glover fiducial u0", lw=4, lc = blues[1])
    plot!(time, C_glo_neu, label="Glover neutral u0", lw=4, lc = blues[4])
    plot!(time, C_umi_fid, label="Umist fiducial u0", lw=4, lc =blues[7])
    plot!(time, C_umi_neu, label="Umist neutral u0", lw=4, lc = blues[9], 
        xlabel="Time (Years)", 
        ylabel="Abundance per Hydrogen Molecule", 
        title="Glover Vs. Umist: Carbon", 
        legend=:topright, 
        xticks=(xtick_vals, xtick_labels), 
        yticks=(ytick_vals, ytick_labels),
        legendfontsize=10, 
        xguidefontsize=14, 
        yguidefontsize=14, 
        ytickfontsize=10,
        xtickfontsize=10,
        fontfamily="Bookman Light", 
        dpi=300)
    savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/scripts/Data/glover_umist_carbon_abundance_plot.svg")
=#

    #=
    ### Glover Vs. Umist: Oxygen ###
    plot(time, O_glo_fid, label="Glover fiducial u0", lw=4, lc = golds[1])
    plot!(time, O_glo_neu, label="Glover neutral u0", lw=4, lc = golds[4])
    plot!(time, O_umi_fid, label="Umist fiducial u0", lw=4, lc = golds[7])
    plot!(time, O_umi_neu, label="Umist neutral u0", lw=4, lc = golds[9], 
          xlabel="Time (Years)", 
          ylabel="Abundance per Hydrogen Molecule", 
          title="Glover Vs. Umist: Oxygen", 
          legend=:bottomright, 
          xticks=(xtick_vals, xtick_labels),
          yticks=(ytick_vals2, ytick_labels2),
          legendfontsize=10, 
          xguidefontsize=14, 
          yguidefontsize=14, 
          ytickfontsize=10,
          xtickfontsize=10,
          fontfamily="Bookman Light", 
          dpi=300)
    savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/scripts/Data/glover_umist_oxygen_abundance_plot.svg")
=#


    ### Relative Error Plot ###
    sheet3 = xf["Sheet3"]
    nel_fid_rel = [sheet3["G$(i)"] for i in 18:31]
    glo_fid_rel = [sheet3["H$(i)"] for i in 18:31]

    nel_neu_rel = [sheet3["O$(i)"] for i in 18:31]
    glo_neu_rel = [sheet3["P$(i)"] for i in 18:31]
    species_names = ["H2", "H3+", "e", "He", "He+", "C", "CH", "CH2", "O", "H2O", "O2", "CO", "HCO+", "C+"]

   
      # Remove the 7th element
      del_species_names = deleteat!(copy(species_names), 7)
      del_nel_fid_rel = deleteat!(copy(nel_fid_rel), 7)
      del_glo_fid_rel = deleteat!(copy(glo_fid_rel), 7)
      scatter(del_species_names, del_nel_fid_rel, label="Nelson fiducial relative error", 
            marker =:circle, 
            markersize = 10, 
            markerstrokewidth = 0.0, 
            markercolor = blues[3])
      #scatter!(species_names, nel_fid_rel, label="Nelson fiducial relative error", marker = (:star5, 6))
      #scatter!(species_names, glo_neu_rel, label="Glover neutral relative error", marker = (:square, 6))
      scatter!(del_species_names, del_glo_fid_rel, label="Glover fiducial relative error", 
            xlabel = "Species",
            ylabel = "Relative Error",
            title = "Relative Error by Species",
            marker =:diamond, 
            legend=:bottomright,
            markersize = 8,
            markerstrokewidth = 0.0,
            yscale = :log10, 
            markercolor = blues[7], 
            fontfamily="Bookman Light",
            xguidefontsize=14, 
            yguidefontsize=14, 
            ytickfontsize=10,
            xtickfontsize=10, 
            dpi=300)
      savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/scripts/Data/rel_error_plot.svg")
  

end
