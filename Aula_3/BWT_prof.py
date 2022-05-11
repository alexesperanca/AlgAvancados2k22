# -*- coding: utf-8 -*-

# Cenas stor a mostrar na aula

seq = "TATAGAGAC$"
perm_ord = sorted([(seq[i:] + seq[:i], i) for i in range(len(seq))])
print(perm_ord)
bwt, suffix_array = zip(*[(s[-1], p) for s, p in perm_ord])
bwt = ''.join(bwt)
print(bwt, suffix_array)

def tabela():                    # Esta função tem outra função imbutida onde tem 1 dicionário dedicado a cada objeto criado por esta função -> útil para obter algo único a cada objeto
	D = {}
	def _add(x):
		nonlocal D               # Permite a utilização do dicionário vazio "D" criado antes
		idx = D.get(x, 0)        # Obtém o valor de "x" do dicionário, mas caso n retorne nada dá "0"
		D[x] = idx + 1
		return x + str(idx)
	return _add

fun = tabela()
bwt_off = [fun(x) for x in bwt]  # Uso do dicionário criado dentro da função e atribuição da contagem obtida à letra
print("Bwt_off:", bwt_off)
fun = tabela()
ord_off = [fun(x) for x in sorted(bwt)]
print("Bwt_ord:", ord_off)

tab = {k:v for k, v in zip(bwt_off, ord_off)}

rec = ""
x = tab["$0"]
while x != "$0":
	rec += x[0]
	x = tab[x]

print(seq, rec)

print([tab[tab[x]] for x in tab if x[0] == "T"])