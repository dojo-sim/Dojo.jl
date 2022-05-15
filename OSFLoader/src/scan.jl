function nerf_scan(nerf_object, xrange, yrange, zrange; title="nerf density")

    anim = @animate for z ∈ zrange
        density = slice_density(nerf_object, xrange, yrange, z)
        plt = heatmap(xrange, yrange, density,
            c=cgrad([:blue, :white,:red, :yellow]),
            xlabel="x", ylabel="y",
            title=title)
        display(plt)
    end

    return anim
end
