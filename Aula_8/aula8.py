import random

kmers = "ATC CCG CGA GAT TCC".split()

seq=''.join(random.choices("ACTG", k=20))

#print(seq)

kmers = [seq[i: i+3] for i in range(len(seq) -2)]

kmers

def get_edges(kmers): 
    for i,u in enumerate(kmers):
        for j, v in enumerate(kmers):
            if i != j and u[1:] == v[:-1]:
                print(f'{u}_{i} -» {v}_{j}')

            

#print(kmers)


#print(get_edges(kmers))


#print(sorted(kmers))

#get_edges(kmers)

P = "ATC TCC CCG CGA GAT".split()
print(P)

print([x[-1] for x in P])

print(''.join([x[-1] for x in P]))

print(P[0][:-1] + ''.join([x[-1] for x in P]))


print('-------------')

print('Transformar em vertices as letras que estão em comum, assim sabemos que podemos passar nos vertices e ir "cortando" os ramos para reconstruir a sequencia original')
for k in kmers:
    print(f'{k[:-1]} -» {k[1:]}')


# com grafos de debruijn é linear, facil de fazer a reconstruçao se tivermos as sequencias certas'

# Fazer dicionario de dicionarios em relacao aos vertices (dicionario "AC" tem dentro um dicionario para as suas saidas (nao entradas) e um valor associado o nº de saidas que tem para essa saida). Ao fazer a reconstrução, vai-se retirar -1 a esse valor e quando for 0 retira-se essa entrada completamente.
