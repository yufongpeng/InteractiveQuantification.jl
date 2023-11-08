module InteractiveCalibration
using GLMakie, GLM, CSV, TypedTables, Plotly
using Gtk4: save_dialog
export Calibration, Project, project, calibration, plot_cal!, inv_predict, view_cal, view_sample

mutable struct Calibration
    type::Bool
    zero::Bool
    weight::Float64
    formula::FormulaTerm
    source::Table
    model
end

mutable struct Project
    calibration::Calibration
    sample::Table
end

include("calibration.jl")
include("calplot.jl")
include("tableview.jl")

end
