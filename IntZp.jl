using LinearAlgebra
import Base: +, -, *, inv, /, isless, iszero
import LinearAlgebra: dot, inv, zeros, ones, *

struct IntZp <: Number
    n::Int
    p::Int
    function IntZp(n, p) 
        a = n % p
    if a < 0
        a += p
    end
        return new(a, p)
    end
end

function +(a::IntZp, b::IntZp)
    a.p != b.p && error("$a e $b não estão no mesmo corpo")
    x = (a.n + b.n) % a.p
    return IntZp(x, a.p)
end

function -(a::IntZp, b::IntZp)
    a.p != b.p && error("$a e $b não estão no mesmo corpo")
    x = (a.n - b.n) % a.p
    return IntZp(x, a.p)
end

function *(a::IntZp, b::IntZp)
    a.p != b.p && error("$a e $b não estão no mesmo corpo")
    x = (a.n * b.n) % a.p
    return IntZp(x, a.p)
end

function inv(num::IntZp)
    x = estendido(num)
    return IntZp(x, num.p)
end

function /(a::IntZp, b::IntZp)
    a.p != b.p && error("$a e $b não estão no mesmo corpo")
    b = inv(b)
    x = a.n * b.n
    return IntZp(x, a.p)
end

function isless(a::IntZp, b::IntZp)
    a.p != b.p && error("$a e $b não estão no mesmo corpo")
    return isless(a.n, b.n)
end

function iszero(a::IntZp)
    return iszero(a.n)
end

function zeros(a::IntZp, n::Int)
    return fill(IntZp(0, a.p), n)
end

function zeros(a::IntZp, n::Int, m::Int)
    return fill(IntZp(0, a.p), n, m)
end

function ones(a::IntZp, n::Int)
    return fill(IntZp(1, a.p), n)
end

function *(n::Int, m::IntZp)
prod = n*m.n
convert = prod%m.p
return IntZp(convert, m.p)
end

*(m::IntZp, n::Int) = n*m

function estendido(num::IntZp)
    a = num.n
    b = num.p
    x = 0
    x1 = 1
    y = 1
    y1 = 0
    r = b
    r1 = a
    while r != 0
        quociente = div(r1, r)
        (r1, r) = (r, r1 - quociente * r)
        (x1, x) = (x, x1 - quociente * x)
        (y1, y) = (y, y1 - quociente * y)
    end
    return x1
end

IntZp(A::AbstractArray{Int}, p::Int) = IntZp.(A, p)

function *(A::Array{IntZp}, B::Array{IntZp})
    p1 = A[1,1].p
    p2 = B[1,1].p
    p1 !=p2 && error("As matrizes não pertencem ao mesmo corpo")
    m,n = size(A)
    r,s = size(B)
    n != r && error("As ordens das matrizes não são compatíveis")
    C = convert(Matrix{Int}, zeros(m,n))
    D = convert(Matrix{Int}, zeros(r,s))
    for i=1:m
        for j = 1:n
            C[i,j] += A[i,j].n
        end
    end
    for i=1:r
        for j = 1:s
            D[i,j] += B[i,j].n
        end
    end 
    AB = C*D
    return IntZp(AB, p1)
end

function getnp(A::AbstractArray{IntZp})
    p = A[1,1].p
    m,n = size(A)
    C = convert(Matrix{Int}, zeros(m,n))
    for i=1:m
        for j = 1:n
            C[i,j] += A[i,j].n
        end
    end
    return C, p
end

getnp(a::IntZp) = a.n, a.p
