get_formula(cal::Calibration) = get_formula(cal.type, cal.zero)
get_formula(type::Bool, zero::Bool) = if type 
    zero ? @formula(y ~ 0 + x) : @formula(y ~ x)
else
    zero ? @formula(y ~ 0 + x + x ^ 2) : @formula(y ~ x + x ^ 2)
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

function calibration(file::String)
    name = basename(file)
    dir = dirname(file)
    config_path = joinpath(dir, "(CONFIG)$name")
    tbl = CSV.read(file, Table)
    config_kw, config_vl = eachline(config_path)
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
    model = lm(f, source[source.include]; wts = source.x[source.include] .^ weight)
    xlevel = unique(source.x)
    source = Table(; id = source.id, level = [findfirst(x -> i == x, xlevel) for i in source.x], y = source.y, x = source.x, x̂ = zeros(Float64, length(source)), accuracy = zeros(Float64, length(source)), include = source.include)
    cal = Calibration(type, zero, weight, f, source, model)
    cal.source.x̂ .= inv_predict(cal, source)
    cal.source.accuracy .= cal.source.x̂ ./ source.x
    cal
end

function refit!(cal::Calibration)
    cal.model = lm(cal.formula, cal.source[cal.source.include]; wts = cal.source.x[cal.source.include] .^ cal.weight)
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