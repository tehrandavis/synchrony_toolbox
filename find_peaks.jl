# taken from Peaks.jl (https://github.com/halleysfifthinc/Peaks.jl)

export maxima,
       minima

for (funcname, comp, notcomp) in ((:maxima, :<, :>),
                         (:minima, :>, :<))

    @eval begin
        function ($funcname)(x::AbstractVector,
                        w::Integer=1,
                        strictbounds::Bool=true)
            w > 0 || throw(ArgumentError("window cannot be negative"))
            xlen = length(x)
            idxs = Int[]

            # There can't be more than one peak every `w` elements
            maxN = strictbounds ? cld(xlen,2*w) : fld(xlen,2*w)
            N = 0

            ww = [ collect(-w:-1); collect(1:w) ]

            sizehint!(idxs,maxN)

            firsti = firstindex(x)
            lasti = lastindex(x)
            i = firsti
            maxI = lasti - w

            if !strictbounds
                while i <= firsti + w
                    xi = x[i]
                    peak = true
                    if ismissing(xi)
                        i = something(findnext(!ismissing, x, i+1), lasti+1)
                    elseif isnan(xi)
                        i = something(findnext(!isnan, x, i+1), lasti+1)
                    else
                        for j in max(i-w, firsti):(i+w)
                            i === j && continue
                            xj = x[j]
                            if coalesce(($comp)(xi, xj), strictbounds) || (strictbounds && (ismissing(xj) || isnan(xj)))
                                peak &= false
                                break # No reason to continue checking if the rest of the elements qualify
                            elseif xi === xj
                                k = findnext(y -> xi !== y, x, j+1)
                                if isnothing(k) # x is constant till the end, not a peak
                                    if !strictbounds
                                        push!(idxs,i)
                                        N += 1
                                        i = lasti+1
                                        peak &= false
                                        break
                                    else
                                        i = lasti+1
                                        peak &= false
                                        break
                                    end
                                else
                                    xk = x[k]
                                    if coalesce(($comp)(xi, xk), strictbounds) || (strictbounds && (ismissing(xk) || isnan(xk))) # x moves towards a peak
                                        peak &= false
                                        break
                                    else # Push new peak here to shift the right number of elements
                                        push!(idxs,i)
                                        N += 1
                                        i = max(k,i+w)
                                        peak &= false
                                        break
                                    end
                                end
                            end
                        end

                        if peak
                            push!(idxs,i)
                            N += 1
                            i += w # There can't be another peak for at least `w` more elements
                        end
                        i += 1
                    end
                end
            end

            i = max(i, firsti + w, something(findfirst(y -> y !== x[firsti], x), firsti))

            while i <= maxI
                xi = x[i]
                peak = true
                if ismissing(xi)
                    i = something(findnext(!ismissing, x, i+1), lasti+1)
                elseif isnan(xi)
                    i = something(findnext(!isnan, x, i+1), lasti+1)
                else
                    for j in ww # For all elements within the window
                        xj = x[i+j]
                        if coalesce(($comp)(xi, xj), strictbounds) || (strictbounds && (ismissing(xj) || isnan(xj)))
                            peak &= false
                            break # No reason to continue checking if the rest of the elements qualify
                        elseif xi === xj
                            xi1 = x[i-1]
                            if coalesce(($notcomp)(xi, xi1), !strictbounds) || (!strictbounds && isnan(xi1))
                                k = findnext(y -> xi !== y, x, i+j+1)
                                if isnothing(k)
                                    if !strictbounds # x is constant till the end, not a peak
                                        push!(idxs,i)
                                        N += 1
                                        i = lasti+1
                                        peak &= false
                                        break
                                    else
                                        i = lasti+1
                                        peak &= false
                                        break
                                    end
                                else
                                    xk = x[k]
                                    if coalesce(($comp)(xi, xk), strictbounds) || (strictbounds && (ismissing(xk) || isnan(xk)))
                                        # x moves towards a peak or ends with a missing or NaN
                                        peak &= false
                                        break
                                    elseif j === 1 # Push first element of plateau as peak here to shift the correct number of elements
                                        push!(idxs,i)
                                        N += 1
                                        i = max(k,i+w)
                                        peak &= false
                                        break
                                    end
                                end
                            else
                                peak &= false
                                break
                            end
                        end
                    end

                    if peak
                        push!(idxs,i)
                        N += 1
                        i += w # There can't be another peak for at least `w` more elements
                    end
                    i += 1
                end
            end

            if !strictbounds
                while i <= lasti
                    xi = x[i]
                    peak = true
                    if ismissing(xi)
                        i = something(findnext(!ismissing, x, i+1), lasti+1)
                    elseif isnan(xi)
                        i = something(findnext(!isnan, x, i+1), lasti+1)
                    else
                        for j in (i-w):min((i+w),lasti)
                            i === j && continue
                            xj = x[j]
                            if coalesce(($comp)(xi, xj), strictbounds) || (strictbounds && (ismissing(xj) || isnan(xj)))
                                peak &= false
                                break # No reason to continue checking if the rest of the elements qualify
                            elseif xi === xj
                                xi1 = x[i-1]
                                if coalesce(($notcomp)(xi, xi1), !strictbounds) || (!strictbounds && isnan(xi1))
                                    k = findnext(y -> xi !== y, x, j+1)
                                    if isnothing(k)
                                        if !strictbounds # x is constant till the end, not a peak
                                            push!(idxs,i)
                                            N += 1
                                            i = lasti+1
                                            peak &= false
                                            break
                                        else
                                            i = lasti+1
                                            peak &= false
                                            break
                                        end
                                    else
                                        xk = x[k]
                                        if coalesce(($comp)(xi, xk), strictbounds) || (strictbounds && (ismissing(xk) || isnan(xk))) # x moves towards a peak
                                            peak &= false
                                            break
                                        else # Push new peak here to shift the right number of elements
                                            push!(idxs,i)
                                            N += 1
                                            i = max(k,i+w)
                                            peak &= false
                                            break
                                        end
                                    end
                                else
                                    peak &= false
                                    break
                                end
                            end
                        end

                        if peak
                            push!(idxs,i)
                            N += 1
                            i += w # There can't be another peak for at least `w` more elements
                        end
                        i += 1
                    end
                end
            end

            resize!(idxs,N)
            sizehint!(idxs,N)

            return idxs
        end
    end
end

@doc """
    maxima(x[, w=1, strictbounds=true])

Find the indices of the local maxima of `x` where each maxima is either the maximum of
`x[-w:w]` or the first index of a plateau.

If `strictbounds` is `true`, all elements of `x[-w:w]` must exist and may not be `missing`
or `NaN`. If `strictbounds` is `false`, elements of `x[-w:w]` may not exist (eg peaks may
be less than `w` indices from either end of `x`), or may be `missing` or `NaN`. `missing`
or `NaN` must not be maxima.
"""
maxima

@doc """
    minima(x[, w=1, strictbounds=false])

Find the indices of the local minima of `x` where each minima is either the minimum of
`x[-w:w]` or the first index of a plateau.

If `strictbounds` is `true`, all elements of `x[-w:w]` must exist and may not be `missing`
or `NaN`. If `strictbounds` is `false`, elements of `x[-w:w]` may not exist (eg peaks may
be less than `w` indices from either end of `x`), or may be `missing` or `NaN`. `missing`
or `NaN` must not be minima.
"""
minima
