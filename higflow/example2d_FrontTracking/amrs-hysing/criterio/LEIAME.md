Diretório de saída do critério de refino.

`ft_escreve_amr_criterio` escreve aqui o `.amr` multinível a cada remalhamento, e
`input/amr.load.domain.yaml` aponta para ele. O arquivo em si é gerado e não é
versionado — mas **o diretório precisa existir**: `fopen` em modo escrita não o
cria, a escrita falha, e o carregador de domínio then dá SEGV em `fscanf` com
FILE* nulo, longe da causa. Foi o que aconteceu ao apagar a pasta numa limpeza.
