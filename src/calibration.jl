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
    elseif cal.zero
        (-β[1] .+ sqrt.(β[1] ^ 2 .+ 4 .* β[2] .* y)) ./ 2 ./ β[2]
    else
        (-β[2] .+ sqrt.(β[2] ^ 2 .+ 4 .* β[3] .* (y .- β[1]))) ./ 2 ./ β[2]
    end
end

function calibration(tbl::Table; 
                    type = true, 
                    zero = false, 
                    weight = 0
                    )
    id = findall(x -> isa(x, Number), tbl.y)
    tbl = tbl[id]
    source = :id in propertynames(tbl) ? tbl : Table((; id = collect(1:length(tbl)), ), tbl)
    f = get_formula(type, zero)
    model = lm(f, source; wts = source.x .^ weight)
    xlevel = unique(source.x)
    source = Table(source; level = [findfirst(x -> i == x, xlevel) for i in source.x], include = trues(length(source)), x̂ = zeros(Float64, length(source)), accuracy = zeros(Float64, length(source)))
    cal = Calibration(type, zero, weight, f, source, model)
    cal.source.x̂ .= inv_predict(cal, source)
    cal.source.accuracy .= cal.source.x̂ ./ source.x
    cal
end

function refit!(cal::Calibration)
    cal.model = lm(cal.formula, cal.source[cal.source.include]; wts = cal.source.x[cal.source.include] .^ cal.weight)
    cal
end

function project(cal::Table, sample::Table; 
        type = true, 
        zero = false, 
        weight = 0)
    cal = calibration(cal; type, zero, weight)
    sample = Table(sample; x̂ = inv_predict(cal, sample))
    Project(cal, sample)
end
