# Standalone hot-path micro-benchmark for ProteinCoLoc refactoring spike 001.
# Uses ONLY stdlib (Statistics.cor) + built-in @allocated/@elapsed. No real images.
# Replicates current src/ implementations faithfully, then compares targeted alternatives.

using Statistics: cor
using Random

Random.seed!(20260628)

# ----------------------------------------------------------------------------
# CURRENT implementations (copied verbatim from src/colocalization.jl)
# ----------------------------------------------------------------------------
function patch_current(img::AbstractMatrix{T}, num_patches::Integer) where T <: Union{Float64, Missing}
    rows, cols = size(img)
    patch_size_x = rows ÷ num_patches
    patch_size_y = cols ÷ num_patches
    trimmed = @view img[1:num_patches*patch_size_x, 1:num_patches*patch_size_y]
    patches = Array{Union{Float64, Missing}, 4}(undef, num_patches, num_patches, patch_size_x, patch_size_y)
    @inbounds for j in 1:num_patches
        y_start = (j-1)*patch_size_y + 1
        for i in 1:num_patches
            x_start = (i-1)*patch_size_x + 1
            for py in 1:patch_size_y, px in 1:patch_size_x
                patches[i, j, px, py] = trimmed[x_start + px - 1, y_start + py - 1]
            end
        end
    end
    return patches
end

function _exclude_zero_current(a::Vector{Union{T,Missing}}, b::Vector{Union{T,Missing}}) where T <: Number
    n = length(a)
    out_a = Vector{T}(undef, n); out_b = Vector{T}(undef, n); k = 0
    @inbounds for i in 1:n
        ai, bi = a[i], b[i]
        (ismissing(ai) || ismissing(bi)) && continue
        (iszero(ai) || isnan(ai) || iszero(bi) || isnan(bi)) && continue
        k += 1; out_a[k] = ai; out_b[k] = bi
    end
    return resize!(out_a, k), resize!(out_b, k)
end

function correlation_current(x::Array{T,4}, y::Array{T,4}; method::Symbol=:pearson) where T <: Union{Float64,Missing}
    cor_dict = Dict(:pearson => cor)          # Dict built every call; returns ::Function (dynamic dispatch)
    cor_func = cor_dict[method]
    num_patches = size(x, 1)
    ρ = zeros(Union{Float64, Missing}, num_patches, num_patches)
    @inbounds for i in 1:num_patches
        for j in 1:num_patches
            a = vec(collect(view(x, i, j, :, :)))   # collect(view) -> alloc, then vec
            b = vec(collect(view(y, i, j, :, :)))
            a, b = _exclude_zero_current(a, b)
            length(a) <= 15 ? ρ[i, j] = missing : ρ[i, j] = cor_func(a, b)
        end
    end
    return ρ
end

# ----------------------------------------------------------------------------
# PROPOSED alternatives (illustrative; quantify the headline wins)
# ----------------------------------------------------------------------------
# A: patch into a dense Float64 4D array (missing carried as NaN sentinel)
function patch_dense(img::AbstractMatrix{Float64}, num_patches::Integer)
    rows, cols = size(img)
    psx = rows ÷ num_patches; psy = cols ÷ num_patches
    patches = Array{Float64,4}(undef, num_patches, num_patches, psx, psy)
    @inbounds for j in 1:num_patches
        y0 = (j-1)*psy
        for i in 1:num_patches
            x0 = (i-1)*psx
            for py in 1:psy, px in 1:psx
                patches[i, j, px, py] = img[x0+px, y0+py]
            end
        end
    end
    return patches
end

# B: correlation with concrete cor, no Dict, no collect(view): single-pass into reused buffers
function correlation_fast(x::Array{Float64,4}, y::Array{Float64,4})
    np = size(x, 1); psx = size(x, 3); psy = size(x, 4)
    cap = psx*psy
    ρ = Matrix{Union{Float64,Missing}}(undef, np, np)
    bufa = Vector{Float64}(undef, cap); bufb = Vector{Float64}(undef, cap)
    @inbounds for i in 1:np, j in 1:np
        k = 0
        for py in 1:psy, px in 1:psx
            ai = x[i,j,px,py]; bi = y[i,j,px,py]
            (iszero(ai) || isnan(ai) || iszero(bi) || isnan(bi)) && continue
            k += 1; bufa[k] = ai; bufb[k] = bi
        end
        ρ[i,j] = k <= 15 ? missing : cor(view(bufa,1:k), view(bufb,1:k))
    end
    return ρ
end

# ----------------------------------------------------------------------------
# Harness
# ----------------------------------------------------------------------------
function timeit(f, args...; reps=20)
    f(args...)                                   # warmup / compile
    GC.gc()
    best = Inf; allocs = 0
    for _ in 1:reps
        t = @elapsed (r = f(args...))
        best = min(best, t)
    end
    allocs = @allocated f(args...)
    return best, allocs
end

function build_image(sz, frac_zero)
    img = rand(Float64, sz, sz)
    nz = round(Int, frac_zero*length(img))
    idx = randperm(length(img))[1:nz]
    img[idx] .= 0.0                              # simulate post-Otsu masking
    return img
end

println("="^70)
println("ProteinCoLoc hot-path micro-benchmark  (Julia $(VERSION))")
println("="^70)

for sz in (512, 1024, 2048)
    for np in (8, 32)
        img_f = build_image(sz, 0.5)
        img_u = convert(Matrix{Union{Float64,Missing}}, img_f)

        # --- patch ---
        tp_cur, ap_cur = timeit(patch_current, img_u, np)
        tp_new, ap_new = timeit(patch_dense, img_f, np)

        xu = patch_current(img_u, np); yu = patch_current(img_u, np)
        xf = patch_dense(img_f, np);   yf = patch_dense(img_f, np)

        # --- correlation ---
        tc_cur, ac_cur = timeit(correlation_current, xu, yu)
        tc_new, ac_new = timeit(correlation_fast, xf, yf)

        println()
        println("IMG $(sz)x$(sz)  num_patches=$np  (patch size $(sz÷np)x$(sz÷np))")
        println("  patch       current: $(round(tp_cur*1e3,digits=3)) ms  $(round(ap_cur/1e6,digits=2)) MB | dense:  $(round(tp_new*1e3,digits=3)) ms  $(round(ap_new/1e6,digits=2)) MB  | speedup $(round(tp_cur/tp_new,digits=2))x  mem $(round(ap_cur/ap_new,digits=2))x")
        println("  correlation current: $(round(tc_cur*1e3,digits=3)) ms  $(round(ac_cur/1e6,digits=2)) MB | fast:   $(round(tc_new*1e3,digits=3)) ms  $(round(ac_new/1e6,digits=2)) MB  | speedup $(round(tc_cur/tc_new,digits=2))x  mem $(round(ac_cur/ac_new,digits=2))x")
    end
end
println()
println("done")
