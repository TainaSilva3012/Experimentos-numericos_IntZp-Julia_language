using Unicode
using  LinearAlgebra
include("IntZp.jl")

## Codificação

function lerarquivo(caminho) 
    pre_string = open(f->read(f, String), caminho)
    string = Unicode.normalize(pre_string) 
    return string
end 

pre_tabela = lerarquivo("precodificacao.txt")
TABELA = split(pre_tabela,"")

function converte_char(letra) 
    posicao = findall(x-> x==letra, TABELA)
    return posicao[1] - 1
end

function converte_texto(texto)
    texto_array = split(texto,"")
    array_convertido = converte_char.(texto_array)
    return array_convertido
end

function concatena_matrix(array::Array)
    tamanho = length(array)
    resto = tamanho % n
    falta = n - resto
    finalzinho = repeat([111],falta) 
    array_completado = [array; finalzinho]
    return reshape(array_completado,n,:)
end

function chave(n::Int,p::Int)
    C= diagm(fill(1, n))
    for j = 1:n
        C[j,j] = rand(1:p-1)
        for i = 1:n
            C[j, i] += rand(0:p-1)
        end
    end
    return IntZp(C, p)
end

## Decodificação

function subdir(L::Matrix{IntZp}, b::Vector{IntZp})
    @assert L[1].p == b[1].p "L e b não estão no mesmo corpo"
    @assert istril(L) "Erro L não é triangular inferior"
    n = size(L, 1)
    y = copy(b)
    for k = 1:n
        for j = 1:k-1
            y[k] -= L[k, j] * y[j]
        end
        y[k] /= L[k, k]
    end
    return y
end

function subrev(U::Matrix{IntZp}, b::Vector{IntZp})
    p = U[1].p
    @assert p == b[1].p "U e b não estão no mesmo corpo"
    @assert istriu(U) "Erro U não é triangular superior"
    n = size(U, 1)
    x = copy(b)
    for k = n:-1:1
        for j = k+1:n
            x[k] -= U[k, j] * x[j]
        end
        if U[k, k] == IntZp(0, p)
            error("U é singular")
        end
        x[k] /= U[k, k]
    end
    return x
end

function LU(A::Matrix{IntZp})
    n = size(A, 1)
    U = copy(A)
    p = A[1, 1].p
    L = IntZp.(diagm(fill(1, n)), p)
    for j = 1:n-1
        U[j,j] == IntZp(0,p) && error("A matriz dada possui pivô 
        				nulo em Z_$p, use PA = LU")	
        for i = j+1:n
            L[i, j] = U[i, j] / U[j, j]
            for k = j:n
                U[i, k] -= L[i, j] * U[j, k]
            end
        end
    end
    return L, U
end

function inv(A::Matrix{IntZp})
    Ainv = copy(A)
    m, n = size(A)
    @assert m == n "A matriz dada não é quadrada"
    L, U = LU(A)
    a = A[1]
    for colA in eachcol(A)
        j = colA.indices[2]
        ej = zeros(a, n)
        ej[j] = IntZp(1, a.p)
        y = subdir(L, ej)
        Ainv[:,j] = subrev(U, y)
    end
    return Ainv
end

function converte_matrix(matrix)
    return reshape(matrix,1,:) 
end

function converte_num(n)
    return TABELA[n+1]
end

function converte_array(array)
    array_texto = converte_num.(array)
    return join(array_texto,"")
end

function escrevearquivo(caminho, texto)
    open(caminho, "w") do io
        write(io, texto)
    end
end
