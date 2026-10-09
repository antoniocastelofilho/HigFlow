# Retomar o FT N=50 com filtro — estado de 09/10/2026 09:40

Escrito antes de o notebook ser fechado para carregar. Tudo que importa está em
disco; nada depende de processo vivo nem desta sessão.

## 1. Onde a corrida está

**FT viscoelástico, λ₁=0,4 β=0,2, N=50, com filtro de ondulações.**
Diretório: `~/HiGFlow-System/HigFlow/higflow/ve/ft-N50-salva`

| | |
|---|---|
| passo | **988 de 2001 (49%)**, t = 1,976 |
| quadro | 80 |
| ponto de salvamento | `output/` — passo 988, marca e `init.yaml` coerentes, sem `.gravando` |
| cópia de segurança | `seguranca/` — passo 988, conjunto completo (20 campos + frente) |
| quadros de superfície | 80, em `VTKS/` |
| perdido no último corte | 1 passo |

O lançador parou o trecho às 08:19 **no piso de 15%, num instante escolhido** —
parada limpa, não queda.

## 2. Por que parou: a bateria não carrega

Medido com a máquina ociosa, três amostras a 10 s:

- `ADP1 online=1` — adaptador conectado
- `charge_now` **travado em 1,299 Ah** (18,7% de 6,955 Ah) desde 08:20
- `current_now` entre **0 e 46 mA** — nada entrando
- `status` alternando **`Full` / `Discharging`** — o campo mentindo, dizendo "Full" a 18,7%
- `voltage_now` **11,182 V**, abaixo do mínimo de projeto (11,360 V)

O controlador decidiu que a bateria está cheia a 18,7%. É condição de hardware.
Fechar o notebook (suspend) pode reiniciar o controlador embutido e resolver —
um desligamento completo tem mais chance que um suspend. Se na volta
`charge_now` ainda estiver travado, o problema não é de carga em si e vale
considerar recalibração do medidor.

**Nunca decidir carga pelo campo `status`.** Use:
```bash
B=/sys/class/power_supply/BAT0
echo "$(cat $B/capacity)% razao $(( $(cat $B/charge_now)*100/$(cat $B/charge_full) ))% $(cat $B/status) $(( $(cat $B/current_now)/1000 )) mA"
```
`capacity` satura em ~80% nesta máquina (referência do driver ≠ `charge_full`);
a **razão** chega a 100%. As portas do lançador já tratam disso.

## 3. O que deve estar vivo, e como ressuscitar cada peça

Conferir:
```bash
ps -eo pid,args | awk '/lanca.sh$|guarda-bateria|ft-fim|orcamento|ondulacao.sh|corrente-ft/ && !/awk/'
```

| peça | o que faz | relançar com |
|---|---|---|
| `ve/ft-N50-salva/lanca.sh` | espera carga, retoma do 988, atravessa cargas | `cd ~/HiGFlow-System/HigFlow/higflow/ve/ft-N50-salva && ./lanca.sh &` |
| `ve/guarda-bateria.sh` | rede de segurança 12/9/5% | `cd ~/HiGFlow-System/HigFlow/higflow/ve && ./guarda-bateria.sh &` |
| `ve/ft-N50-salva/orcamento-filtro.sh` | custo acumulado do filtro | idem, no diretório da corrida |
| `ve/ft-N50-salva/vigia-ondulacao.sh` | ondulação retida, critério absoluto 4% | idem |
| `ve/corrente-ft-varredura.sh` | dispara a varredura quando o FT fechar | `cd ~/HiGFlow-System/HigFlow/higflow/ve && ./corrente-ft-varredura.sh &` |

O lançador é o único **necessário**: os outros observam. Se só ele subir, a
corrida termina igual.

Registros: `/tmp/ft-N50-status`, `/tmp/ft-N50-salva.log`,
`/tmp/ondulacao-ft.log`, `/tmp/orcamento-filtro.log`, `/tmp/guarda-bateria.log`.
**Estão em /tmp e somem num reboot.** Os de `/tmp/ft-N50-salva.log` em diante
não são necessários para retomar — só para acompanhar.

## 4. Cronograma de retorno

Com a carga cheia e nada competindo, no ritmo medido de **108 s/passo**:

| etapa | passos | CPU | observação |
|---|---|---|---|
| retomar e fechar o FT | 988 → 2001 (1013) | **~30 h** | mais ~35 min de partida por retomada |
| 3º quadro de campo | frame 160, t=4 | ~20 min | 310 MB |
| varredura de cadência | 4 × 2001 em N=25 | **14–17 h** | dispara sozinha depois do FT |

Dreno medido neste caso: **~1 ponto de razão por hora**. De 100% até o piso de
15% dá ~12 h de corrida por carga, logo **2 a 3 paradas** até fechar o FT.
Relógio estimado: **2 a 3 dias** para o FT, mais 1 dia para a varredura.

## 5. Como verificar que retomou certo

```bash
cd ~/HiGFlow-System/HigFlow/higflow/ve/ft-N50-salva
tail -5 /tmp/ft-N50-status                      # deve dizer "RETOMANDO do passo 988, faltam 1013"
LC_ALL=C grep -a 'Reloading properties' /tmp/ft-N50-salva.log | tail -1
LC_ALL=C grep -a 'retomada: frente lida' /tmp/ft-N50-salva.log | tail -1   # nv, nt, chi=2
LC_ALL=C grep -a 'should have been saved' /tmp/ft-N50-salva.log            # VAZIO = certo
ls VTKS/*frente_* | wc -l                        # deve continuar de 80, sem salto
```
A linha `should have been saved` é a que denunciava a retomada quebrada; se ela
aparecer, o `init.yaml` não está sendo gravado e a próxima parada perde o trecho.

## 6. O que já está medido (não refazer)

- **Ondulação**: 0,37% das arestas acima do limiar no quadro 79, máximo 0,055.
  O gêmeo **sem** filtro, na mesma janela até t=1,6: 8,20%, razão **25×**.
  As corridas sem filtro chegaram a 21% com máximo 0,997.
- **A cúspide se formou e o filtro não a tocou**: assimetria cruzou 1 em t≈1,42
  e está em 1,118; R/Req caiu de 0,92 a 0,45. Em t=1,4 filtrado e gêmeo
  concordam em **0,05%** na assimetria e ~1% no R/Req.
- **Integridade**: 80 superfícies fechadas, χ=2 em todos os quadros e nas 49
  aplicações do filtro; volume das aplicações soma 6,6e-15.
- **Deriva de volume da corrida**: 0,71% em t=1,98 (o VOF deu 0,45% em t=4).
- O VOF N=50 **fechou** em 07/10: 2001 passos, z_c(4)=2,0417, 161 quadros.

## 7. Decisões pendentes (suas)

1. **A diferença FT‑VOF de esfericidade (0,0023) precisa ser remedida** com as
   corridas filtradas. A medida sem filtro estava deprimida pela área espúria do
   zigue-zague — o filtro repara essa medida, não a contamina, e portanto o
   0,0023 publicado inclui erro de malha.
2. **Cadência do filtro**: a varredura 0/20/10/5 responde quantitativamente. O
   critério de decisão está escrito em `varredura-cadencia/compara.sh` **antes**
   do resultado.
3. **Se a bateria não voltar a carregar**: baixar o piso do lançador daria
   ~13 h de corrida nos 18% disponíveis, com o guarda parando aos 5%. Não fiz
   por conta própria — envolve o hardware.
