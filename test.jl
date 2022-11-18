function f(x,y)
    if x > abs(y)
        return result = 5 * (9*((x)^2) + 16((y)^2))^(0.5)
    end
    return result = 9*x + 16*abs(y)
end

# min = 20
# mini = 0
# for i in -400:400
#     global min
#     global mini
#     i /= 200
#     fx = f(112/225-i, (-7/25)+i)
#     if fx < min
#         min = fx
#         mini = i
#     end
# end
# println(min)
# println(mini)
# println((112/225-mini) > abs((-7/25+mini)))



∇f = x -> [225*x[1] ; 400*x[2]] .* (1/f(x[1],x[2]))
∇f([16/9 , 1])

function minimise(x)
    mint = 0
    min = 1000
    delta = ∇f(x)
    for t in 1:2000000
        t /= 50000
        output = x - t .* delta
        output = f(output[1],output[2])
        if output < min
            min = output
            mint = t
        end
    end
    return x - mint .* delta
end

x = minimise(minimise([16/9 ; 1]))
x[1]/x[2]
# ∇f([112/225; -7/25])
