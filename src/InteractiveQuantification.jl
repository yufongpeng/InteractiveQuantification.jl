module InteractiveQuantification
using GLMakie, GLM, CSV, TypedTables, Plotly, LinearAlgebra, Reexport
using Gtk4: save_dialog
export plot!, view_cal, view_sample
@reexport using ChemistryQuantitativeAnalysis

include("calplot.jl")
include("tableview.jl")

end
