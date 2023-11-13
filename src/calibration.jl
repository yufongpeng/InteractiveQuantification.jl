get_formula(cal::Calibration) = get_formula(cal.type, cal.zero)
get_formula(type::Bool, zero::Bool) = if type 
    zero ? @formula(y ~ 0 + x) : @formula(y ~ x)
else
    zero ? @formula(y ~ 0 + x + x ^ 2) : @formula(y ~ x + x ^ 2)
end

inv_predict(project::Project, tbl) = inv_predict(project.calibration, tbl)
function inv_predict_sample!(project::Project)
    project.sample.x̂ .= inv_predict(project.calibration, project.sample)
    project
end
inv_predict_cal!(project::Project) = (inv_predict_cal!(project.calibration); project)
function inv_predict_cal!(cal::Calibration)
    cal.source.x̂ .= inv_predict(cal, cal.source)
    cal
end
inv_predict(cal::Calibration, tbl::Table) = inv_predict(cal, getproperty(tbl, cal.formula.lhs.sym))
function inv_predict(cal::Calibration, y::AbstractArray)
    β = cal.model.model.pp.beta0
    if cal.type && cal.zero
        y ./ β[1]
    elseif cal.type
        (y .- β[1]) ./ β[2]
    else
        c, b, a = cal.zero ? (0, β...) : β
        d = @. max(b ^ 2 + 4 * a * (y - c), 0)
        @. (-b + sqrt(d)) / 2a
    end
end

accuracy(project::Project, tbl = project.calibration.source) = accuracy(project.calibration, tbl)
accuracy(cal::Calibration, tbl = cal.source) = accuracy(inv_predict(cal, tbl), tbl.x)
accuracy!(project::Project) = (accuracy!(project.calibration); project)
function accuracy!(cal::Calibration)
    cal.source.accuracy .= accuracy(cal.source.x̂, cal.source.x)
    cal
end
accuracy(x̂::AbstractVector, x::AbstractVector) = @. x̂ / x

inv_predict_accuracy! = accuracy! ∘ inv_predict_cal!

function calibration(file::String)
    tbl = CSV.read(joinpath(file, "calibration.csv"), Table)
    config_kw, config_vl = eachline(joinpath(file, "config.csv"))
    config = Dict{Symbol, Any}()
    for (k, v) in zip(split(config_kw, ","), split(config_vl, ","))
        if k == "type" || k == "zero"
            v = v == "TRUE" || v == "True" || v == "true" || v == ""
        elseif k == "weight"
            v = parse(Float64, v)
        else
            continue
        end
        get!(config, Symbol(k), v)
    end
    calibration(tbl; config...)
end

function calibration(tbl::Table; 
                    type = true, 
                    zero = false, 
                    weight = 0
                    )
    id = findall(x -> isa(x, Number), tbl.y)
    tbl = tbl[id]
    source = :id in propertynames(tbl) ? tbl : Table((; id = collect(1:length(tbl)), ), tbl)
    source = :include in propertynames(tbl) ? source : Table(source; include = trues(length(source)))
    f = get_formula(type, zero)
    model = calfit(source, f, type, zero, weight)
    xlevel = unique(source.x)
    source = Table(; id = source.id, level = [findfirst(x -> i == x, xlevel) for i in source.x], y = source.y, x = source.x, x̂ = zeros(Float64, length(source)), accuracy = zeros(Float64, length(source)), include = source.include)
    cal = Calibration(type, zero, weight, f, source, model)
    inv_predict_accuracy!(cal)
    cal
end

function calfit(tbl::Table, formula, type, zero, weight)
    model = lm(formula, tbl[tbl.include]; wts = tbl.x[tbl.include] .^ weight)
    if !type && !zero && model.model.pp.beta0[1] == 0
        m = hcat(ones(Float64, count(tbl.include)), tbl.x[tbl.include], tbl.x[tbl.include] .^ 2)
        sqrtw = diagm(sqrt.(tbl.x[tbl.include] .^ weight))
        y = tbl.y[tbl.include]
        model.model.pp.beta0 = (sqrtw * m) \ (sqrtw * y)
    end
    model
end

function calfit!(cal::Calibration)
    cal.model = calfit(cal.source, cal.formula, cal.type, cal.zero, cal.weight)
    cal
end

project(cal::Table, sample = ""; 
        type = true, 
        zero = false, 
        weight = 0) = project(calibration(cal; type, zero, weight), sample)

project(cal::String, sample = "") = project(calibration(cal), sample)

function project(cal::Calibration, sample = "")
    sample = if sample == ""
        tbl = Table(; id = String[], y = Float64[])
        Table(tbl; x̂ = inv_predict(cal, tbl))
    elseif sample isa String
        tbl = CSV.read(sample, Table)
        Table(tbl; x̂ = inv_predict(cal, tbl))
    else
        sample
    end
    Project(cal, sample)
end

function cretical_point(cal::Calibration)
    β = cal.model.model.pp.beta0
    c, b, a = cal.zero ? (0, β...) : β
    -b / 2a
end

cal_range(project::Project) = cal_range(project.calibration)
cal_range(cal::Calibration) = (lloq(cal), hloq(cal))
lloq(project::Project) = lloq(project.calibration)
hloq(project::Project) = hloq(project.calibration)
lloq(cal::Calibration) = (cal.type || last(cal.model.model.pp.beta0) < 0) ? cal.source.x[findfirst(cal.source.include)] : max(cal.source.x[findfirst(cal.source.include)], cretical_point(cal))
hloq(cal::Calibration) = (cal.type || last(cal.model.model.pp.beta0) > 0) ? cal.source.x[findlast(cal.source.include)] : min(cal.source.x[findlast(cal.source.include)], cretical_point(cal))