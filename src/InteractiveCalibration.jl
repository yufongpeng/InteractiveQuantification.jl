module InteractiveCalibration
using GLMakie, GLM, CSV, TypedTables, Plotly, LinearAlgebra
using Gtk4: save_dialog
export Calibration, Project, project, calibration, plot_cal!, 
    cal_range, lloq, hloq, accuracy, accuracy!,
    inv_predict, inv_predict_cal!, inv_predict_sample!, inv_predict_accuracy!,
    view_cal, view_sample,
    formula_repr, weight_repr, formula_repr_utf8, weight_repr_utf8
    
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
