function plot_cal!(project::Project;
                    lloq_multiplier = 4//3, dev_acc = 0.15,
                    fig_attr = Dict{Symbol, Any}(), 
                    axis_attr = Dict(:title => "Analyte", :xlabel => "Concentration (nM)", :ylabel => "Abundance", :titlesize => 20), 
                    plot_attr = Dict(:point => Dict(:color => [:blue, :red], :inspector_label => (self, i, p) -> string("id: ", project.calibration.source.id[i], "\naccuracy: ", round(project.calibration.source.accuracy[i]; sigdigits = 4))), :line => Dict(:color => :chartreuse)))
    fig = Figure(; fig_attr...)
    menu_type = Menu(fig, options = ["linear", "quadratic"], default = project.calibration.type ? "linear" : "quadratic")
    menu_zero = Menu(fig, options = ["ignore (0, 0)", "include (0, 0)"], default = project.calibration.zero ? "include (0, 0)" : "ignore (0, 0)")
    default_w = weight_repr(project.calibration.weight)
    menu_wt = Menu(fig, options = ["none", "1/√x", "1/x", "1/x²"], default = default_w)
    menu_point = Menu(fig, options = string.(0:length(unique(project.calibration.source.x))), default = "0")
    label_r2 = Label(fig, "R² = $(round(r2(project.calibration.model); sigdigits = 4))"; halign = :left)
    label_formula = Label(fig, formula_repr(project.calibration); halign = :left)
    menu_obj = Menu(fig, options = ["Cal", "Sample", "Fig"], default = "Cal"; width = 70)
    button_show = Button(fig, label = "show")
    button_save = Button(fig, label = "save")
    fig[1, 2] = vgrid!(
        label_r2,
        label_formula,
        Label(fig, "Zoom", width = nothing),
        menu_point,
        Label(fig, "Type", width = nothing),
        menu_type, 
        Label(fig, "Zero", width = nothing),
        menu_zero,
        Label(fig, "Weight", width = nothing),
        menu_wt,
        hgrid!(menu_obj, button_show, button_save);
        tellheight = false
    )
    ax = Axis(fig[1, 1]; axis_attr...)
    scs = scatter!(ax, project.calibration.source.x, project.calibration.source.y; get_point_attr(plot_attr, project.calibration.source.include)...)
    DataInspector(scs)
    #=
    for r in project.calibration.source
        sc = scatter!((r.x, r.y); get_point_attr(plot_attr, Val{r.include})...)
        DataInspector(sc)
        push!(scs, sc)
    end
    =#
    xlevel = unique(project.calibration.source.x)
    xscale = -reduce(-, extrema(xlevel))
    yscale = -reduce(-, extrema(project.calibration.source.y))
    xrange = Table(; x = LinRange(extrema(xlevel)..., convert(Int, reduce(-, extrema(xlevel)) ÷ maximum(xlevel[1:end - 1] .- xlevel[2:end]) * 100)))
    ln = lines!(ax, xrange.x, predict(project.calibration.model, xrange); get!(plot_attr, :line, Dict(:color => :chartreuse))...)
    display(view_cal(project.calibration.source; lloq_multiplier, dev_acc))
    #display(view_sample(project.sample; lloq = project.calibration.source.x[findfirst(project.calibration.source.include)], hloq = project.calibration.source.x[findlast(project.calibration.source.include)], lloq_multiplier, dev_acc))
    # Main.vscodedisplay(project.calibration.source[project.calibration.source.include])
    # fig[1, 3] = vgrid!(map(s -> Label(fig, s; halign = :left), split(sprint(showtable, project.calibration.source), "\n"))...; tellheight = false, width = 250)
    on(events(ax).mousebutton) do event
        if event.action == Mouse.press
            plot, id = pick(ax)
            if id != 0 && plot == scs
                if event.button == Mouse.left
                    project.calibration.source.include[id] = !project.calibration.source.include[id]
                    for (k, v) in pairs(get_point_attr(plot_attr, project.calibration.source.include))
                        getproperty(scs, k)[] = v
                    end
                    refit!(project.calibration)
                    ln.input_args[2][] = predict(project.calibration.model, xrange)
                    project.calibration.source.x̂ .= inv_predict(project.calibration, project.calibration.source)
                    project.calibration.source.accuracy .=  project.calibration.source.x̂ ./ project.calibration.source.x
                    label_r2.text = "R² = $(round(r2(project.calibration.model); sigdigits = 4))"
                    label_formula.text = formula_repr(project.calibration)
                    project.sample.x̂ .= inv_predict(project.calibration, project.sample)  
                    # update_leaf!(tbl, project.calibration.source)
                end
            end
        end
        return Consume(false)
    end
    on(menu_type.selection) do s
        project.calibration.type = s == "linear"
        project.calibration.formula = get_formula(project.calibration)
        refit!(project.calibration)
        ln.input_args[2][] = predict(project.calibration.model, xrange)
        project.calibration.source.x̂ .= inv_predict(project.calibration, project.calibration.source)
        project.calibration.source.accuracy .=  project.calibration.source.x̂ ./ project.calibration.source.x
        label_r2.text = "R² = $(round(r2(project.calibration.model); sigdigits = 4))"
        label_formula.text = formula_repr(project.calibration)
        project.sample.x̂ .= inv_predict(project.calibration, project.sample)
        # update_leaf!(tbl, project.calibration.source)
    end
    on(menu_zero.selection) do s
        project.calibration.zero = s == "include (0, 0)"
        project.calibration.formula = get_formula(project.calibration)
        refit!(project.calibration)
        ln.input_args[2][] = predict(project.calibration.model, xrange)
        project.calibration.source.x̂ .= inv_predict(project.calibration, project.calibration.source)
        project.calibration.source.accuracy .= project.calibration.source.x̂ ./ project.calibration.source.x
        label_r2.text = "R² = $(round(r2(project.calibration.model); sigdigits = 4))"
        label_formula.text = formula_repr(project.calibration)
        project.sample.x̂ .= inv_predict(project.calibration, project.sample)
        # update_leaf!(tbl, project.calibration.source)
    end
    on(menu_wt.selection) do s
        project.calibration.weight = weight_value(s)
        refit!(project.calibration)
        ln.input_args[2][] = predict(project.calibration.model, xrange)
        project.calibration.source.x̂ .= inv_predict(project.calibration, project.calibration.source)
        project.calibration.source.accuracy .=  project.calibration.source.x̂ ./ project.calibration.source.x
        label_r2.text = "R² = $(round(r2(project.calibration.model); sigdigits = 4))"
        label_formula.text = formula_repr(project.calibration)
        project.calibration.source.x̂ .= inv_predict(project.calibration, project.calibration.source)
        project.calibration.source.accuracy .=  project.calibration.source.x̂ ./ project.calibration.source.x
        project.sample.x̂ .= inv_predict(project.calibration, project.sample)
        # update_leaf!(tbl, project.calibration.source)
    end
    on(menu_point.selection) do s
        s = parse(Int, s)
        if s == 0
            autolimits!(ax)
        else
            x_value = xlevel[s] 
            id = findall(==(x_value), project.calibration.source.x)
            y_value = project.calibration.source.y[id]
            Δy = -reduce(-, extrema(y_value))
            yl = extrema(y_value) .+ (-Δy, Δy)
            Δx = Δy * xscale / yscale
            xl = x_value .+ (-Δx, Δx)
            limits!(ax, xl, yl)
        end
    end
    on(button_show.clicks) do s
        if menu_obj.selection[] == "Fig"
            return
        elseif menu_obj.selection[] == "Cal"
            display(view_cal(project.calibration.source; lloq_multiplier, dev_acc))
        else
            display(view_sample(project.sample; lloq = project.calibration.source.x[findfirst(project.calibration.source.include)], hloq = project.calibration.source.x[findlast(project.calibration.source.include)], lloq_multiplier, dev_acc))
        end
        #Main.vscodedisplay(project.calibration.source[project.calibration.source.include])
    end
    on(button_save.clicks) do s
        if menu_obj.selection[] == "Fig"
            save_dialog("Save as", nothing, ["*.png"]; start_folder = pwd()) do f
                save(f, fig; update = false)
            end
            return
        end
        save_dialog("Save as", nothing, ["*.csv"]; start_folder = pwd()) do f
            if menu_obj.selection[] == "Cal"
                CSV.write(f, project.calibration.source)
            else
                CSV.write(f, project.sample)
            end
        end
    end
    fig
end

get_point_attr(plot_attr::Dict, incl::Bool) = NamedTuple(k => incl ? v[1] : v[2] for (k, v) in get!(plot_attr, :point, Dict(:color => [:blue, :red])))
get_point_attr(plot_attr::Dict, incl::BitVector) = NamedTuple(k => isa(v, Vector) ? map(inc -> inc ? v[1] : v[2], incl) : v for (k, v) in get!(plot_attr, :point, Dict(:color => [:blue, :red])))

function weight_repr(cal::Calibration)
    cal.weight in [-0.5, -1, -2] || (cal.weight = 0)
    weight_repr(cal.weight)
end
weight_repr(weight::Number) = if weight == -0.5
    "1/√x"
elseif weight == -1
    "1/x"
elseif weight == -2
    "1/x²"
else
    "none"
end

weight_value(weight) = if weight == "1/√x"
    -0.5
elseif weight == "1/x"
    -1
elseif weight == "1/x²"
    -2
else
    0
end

function formula_repr(cal::Calibration)
    β = cal.model.model.pp.beta0
    if cal.type 
        cal.zero ? "y = $(round(β[1]; sigdigits = 4))x" : "y = $(round(β[1]; sigdigits = 4)) + $(round(β[1]; sigdigits = 4))x"
    else
        cal.zero ? "y = $(round(β[1]; sigdigits = 4))x + $(round(β[1]; sigdigits = 4))x²" : "y = $(round(β[1]; sigdigits = 4)) + $(round(β[2]; sigdigits = 4))x + $(round(β[3]; sigdigits = 4))x²"
    end
end
