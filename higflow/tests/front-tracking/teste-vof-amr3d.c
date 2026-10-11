// Oraculos do criterio de refino do VOF 3D -- a parte que nao precisa do solver.
//
// O QUE SE TESTA AQUI, e por que estas e nao outras:
//
//  1. O PISO DE CINCO VENCE A BANDA.  O PLIC reconstroi um plano por celula a
//     partir da fracao na vizinhanca, e o transporte da fracao precisa da mesma
//     vizinhanca: com menos de cinco celulas finas por lado a reconstrucao nao
//     fecha e o transporte perde massa.  Obedecer em silencio a uma banda que
//     pede menos de cinco quebraria o metodo, e o sintoma -- perda de massa --
//     apareceria longe da causa.  Este e' o oraculo que o usuario pediu.
//
//  2. A BANDA MANDA quando pede cinco ou mais, e a contagem de celulas CRESCE
//     ao refinar.  E' o sentido certo: largura fisica constante entre
//     resolucoes, contagem variavel -- o contrario do que `celmin` fixo faz.
//
//  3. SEM a variavel, nada muda.  Um padrao que silenciosamente alterasse
//     `celmin` mudaria toda corrida existente sem aviso.
//
//  4. Banda invalida (zero, negativa, texto) NAO altera `celmin`.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef double real;

// A funcao sob teste, na forma em que vive em vof-amr-3d.c.  O print0f de la'
// depende do solver; aqui ele e' um printf, e e' a unica diferenca.
static int celmin_da_banda(const char *var, int celmin, real lx, int nx, int niveis)
{
    const char *b = getenv(var);
    if (b == NULL) return celmin;
    real banda = atof(b);
    if (!(banda > 0.0)) return celmin;
    real h_fino = lx / (real) nx;
    for (int l = 0; l < niveis; l++) h_fino *= 0.5;
    int novo = (int) (banda / h_fino + 0.5);
    if (novo < 5) {
        printf("        (aviso) %s=%g pede %d celulas em h=%.5f; o piso de 5 vence\n",
               var, (double) banda, novo, (double) h_fino);
        novo = 5;
    }
    return novo;
}

static int falhas = 0;
static void ok(int c, const char *msg, const char *det)
{
    printf("  [%s] %-50s %s\n", c ? "ok" : "FALHA", msg, det);
    if (!c) falhas++;
}

int main(void)
{
    char d[128];
    printf("\n=== criterio de refino do VOF 3D: a banda e o piso de cinco ===\n\n");

    // --- 1. o piso de cinco vence uma banda fina demais
    setenv("T_BANDA", "0.01", 1);              // h_fino = 1/50 = 0,02 -> pediria 0 ou 1
    int c = celmin_da_banda("T_BANDA", 7, 1.0, 25, 1);
    snprintf(d, sizeof d, "banda 0,01 em h_fino 0,02 -> celmin=%d", c);
    ok(c == 5, "o PISO DE CINCO vence banda fina demais", d);

    setenv("T_BANDA", "0.06", 1);              // pediria 3
    c = celmin_da_banda("T_BANDA", 7, 1.0, 25, 1);
    snprintf(d, sizeof d, "banda 0,06 -> pediria 3, saiu %d", c);
    ok(c == 5, "o piso vence tambem quando a banda pede 3", d);

    // --- 2. a banda manda quando pede cinco ou mais
    setenv("T_BANDA", "0.10", 1);              // 0,10/0,02 = 5 exato
    c = celmin_da_banda("T_BANDA", 7, 1.0, 25, 1);
    snprintf(d, sizeof d, "banda 0,10 em h_fino 0,02 -> celmin=%d", c);
    ok(c == 5, "banda que pede exatamente cinco da' cinco", d);

    setenv("T_BANDA", "0.30", 1);
    c = celmin_da_banda("T_BANDA", 5, 1.0, 25, 1);
    snprintf(d, sizeof d, "banda 0,30 em h_fino 0,02 -> celmin=%d", c);
    ok(c == 15, "banda larga manda sobre o celmin dado", d);

    // --- 3. A MESMA banda em duas resolucoes: a contagem CRESCE ao refinar.
    //        E' o ponto todo da banda fisica -- com celmin fixo, a largura
    //        fisica cairia a' metade e a esteira ficaria de fora.
    setenv("T_BANDA", "0.20", 1);
    const int c_grosso = celmin_da_banda("T_BANDA", 5, 1.0, 25, 1);   // h_fino 0,02
    const int c_fino   = celmin_da_banda("T_BANDA", 5, 1.0, 50, 1);   // h_fino 0,01
    snprintf(d, sizeof d, "h_fino 0,02 -> %d celulas;  h_fino 0,01 -> %d",
             c_grosso, c_fino);
    ok(c_fino == 2 * c_grosso, "mesma banda fisica: a contagem DOBRA ao refinar", d);

    // --- 4. sem a variavel, nada muda
    unsetenv("T_BANDA");
    c = celmin_da_banda("T_BANDA", 7, 1.0, 25, 1);
    snprintf(d, sizeof d, "celmin dado 7, saiu %d", c);
    ok(c == 7, "sem a variavel, celmin passa intacto", d);

    // --- 5. banda invalida nao altera celmin
    int mal = 0;
    const char *ruins[] = {"0", "-1", "abc", ""};
    for (int i = 0; i < 4; i++) {
        setenv("T_BANDA", ruins[i], 1);
        if (celmin_da_banda("T_BANDA", 7, 1.0, 25, 1) != 7) mal++;
    }
    snprintf(d, sizeof d, "quatro entradas invalidas, %d alteraram celmin", mal);
    ok(mal == 0, "banda invalida nao altera celmin", d);

    // --- 6. GUARDA DE VACUIDADE: o teste acima so' vale se a funcao de fato
    //        MUDA celmin em algum caso.  Uma funcao que devolvesse sempre o
    //        argumento passaria em 3, 4, 5 e 6.
    setenv("T_BANDA", "0.40", 1);
    const int mudou = celmin_da_banda("T_BANDA", 5, 1.0, 25, 1);
    snprintf(d, sizeof d, "celmin 5 -> %d", mudou);
    ok(mudou != 5, "a funcao de fato altera celmin em algum caso", d);

    printf("\n%s  (%d falha(s))\n\n", falhas == 0 ? "PASSOU" : "FALHOU", falhas);
    return falhas != 0;
}
