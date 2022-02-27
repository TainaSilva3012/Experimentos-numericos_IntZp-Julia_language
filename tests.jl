using Random

include("auxiliar_functions.jl")
Random.seed!(23)
n = 7
p = 113

##Codificação

chave_codifica = chave(n,p)
poema = lerarquivo("poema.txt")
vetor_palavras = converte_texto(poema)
matriz_concatenada = concatena_matrix(vetor_palavras)
matriz_concatenada = IntZp.(matriz_concatenada, 113)
texto_codificado = chave_codifica * matriz_concatenada
array_codificado = converte_matrix(texto_codificado)
array_codificado_int,  = getnp(array_codificado)
escrevearquivo("poema_codificado.txt",converte_array(array_codificado_int))

## Decodificação

chave_decodifica = inv(chave_codifica)
matriz_concatenada_rev = chave_decodifica * texto_codificado
array_decodificado = converte_matrix(matriz_concatenada_rev)
array_int, p = getnp(array_decodificado)
escrevearquivo("poema_decodificado.txt",converte_array(array_int))
