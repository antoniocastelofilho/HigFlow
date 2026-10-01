// ---------------------------------------------------------------------------
// Fracao volumetrica INICIAL em 3D, por subdivisao subcelular.
//
// POR QUE ESTE ARQUIVO EXISTE.
//
// Os exemplos 2D nao amostram a forma inicial no centro da celula: eles
// SUBDIVIDEM a celula de interface em 16x16 subcelulas, resolvem cada uma
// analiticamente e condensam a soma na celula grossa (ver
// example2d_VOF/ns-example-2d.c, laco do get_fracvolN).  O example3d_VOF_rising
// nascera' com uma rampa linear de uma celula no lugar disso, e o erro no
// volume da bolha, MEDIDO contra 4/3 pi R^3, era:
//
//     N=16   rampa 1.23e-02   subdivisao 16^3  6.11e-05     (200x menor)
//     N=32         4.09e-03                    1.53e-05
//     N=64         5.48e-04                    3.82e-06
//     N=128        2.76e-04                    9.54e-07
//
// A rampa nao e' "de primeira ordem": a ordem medida pula 1.59, 2.90, 0.99.
// Ela nao converge, porque a largura da rampa esta' amarrada a h e o padrao
// de amostragem contra a esfera muda a cada resolucao.  E 1.2% de massa
// errada em N=16 fica errada para sempre: o VOF conserva o que recebeu.
//
// COMO A FRACAO DE CADA SUBCELULA SAI.
//
// Amostra-se `func3` nos oito cantos da subcelula.  Todos do mesmo sinal, a
// subcelula esta' cheia ou vazia.  Senao, ajusta-se o PLANO que a parte
// linear do interpolante trilinear define -- valor medio dos oito cantos e
// gradiente das diferencas -- e aplica-se a forma fechada de Scardovelli &
// Zaleski para o volume de uma caixa cortada por um plano.
//
// Consequencia que serve de oraculo afiado: se `func3` for LINEAR (um plano),
// o ajuste e' exato e a fracao sai correta na precisao da maquina, em UMA
// celula, sem subdividir.  O arreio exige isso.
// ---------------------------------------------------------------------------

// Incluido por #include, nao compilado a' parte -- mesma convencao do
// vof-geometry-2d.c: os Makefiles dos exemplos so' listam ns-example-3d.c, e os
// tipos (`real`, `Point`) vem dos includes que o exemplo ja' fez.  O arreio
// standalone em tests/vof-3d/ declara esses tipos por conta propria.
#ifndef HIGFLOW_EXAMPLES_VOF_GEOMETRY_3D
#define HIGFLOW_EXAMPLES_VOF_GEOMETRY_3D

#include <math.h>

// Level-set da forma inicial: > 0 DENTRO da fase 1 (fracvol = 1).  Mesma
// convencao de sinal do `func` dos exemplos 2D.  Cada exemplo define a sua.
real func3(Point p);

// VOLUME/AREA/COMPRIMENTO do lado NEGATIVO de um plano, em 3, 2 e 1 eixos.
// As tres tem a mesma inclusao-exclusao; existem separadas porque o
// denominador da 3D e' 6 n0 n1 n2 e ele VAI A ZERO quando a normal se alinha
// a um eixo.  Nao e' detalhe de precisao: com duas componentes em 1e-9 o
// denominador cai a 6e-18 enquanto os termos do numerador ficam em 1e-6 --
// dezesseis digitos de cancelamento, nada sobra.  Aparar a normal nao cura,
// causa.  Reduzir a dimensao cura: normal com uma componente nula corta um
// prisma, com duas, uma fatia.  Exige n_i > 0 nos eixos que ficam.
static real _sz3(const real n[3], real alfa, const real L[3]) {
    real den = 6.0 * n[0] * n[1] * n[2];
    real soma = alfa*alfa*alfa;
    for (int i = 0; i < 3; i++) {
        real m = alfa - n[i]*L[i];
        if (m > 0.0) soma -= m*m*m;
    }
    for (int i = 0; i < 3; i++)
        for (int j = i+1; j < 3; j++) {
            real w = alfa - n[i]*L[i] - n[j]*L[j];
            if (w > 0.0) soma += w*w*w;
        }
    real w = alfa - n[0]*L[0] - n[1]*L[1] - n[2]*L[2];
    if (w > 0.0) soma -= w*w*w;
    return soma / den;
}
static real _sz2(real n0, real n1, real alfa, real L0, real L1) {
    real soma = alfa*alfa, m;
    m = alfa - n0*L0;         if (m > 0.0) soma -= m*m;
    m = alfa - n1*L1;         if (m > 0.0) soma -= m*m;
    m = alfa - n0*L0 - n1*L1; if (m > 0.0) soma += m*m;
    return soma / (2.0*n0*n1);
}
static real _sz1(real n0, real alfa, real L0) {
    real x = alfa/n0;
    return (x < 0.0) ? 0.0 : ((x > L0) ? L0 : x);
}

// Fracao de UMA caixa.  Os oito valores de canto dao a DIRECAO do plano; o
// valor no centro (`fcen`) da a sua POSICAO -- ver o comentario sobre a
// curvatura adiante, que e' o que separa 6.1e-05 de 1.2e-04.
// `f[8]` na ordem (i,j,k) com k variando mais rapido:
//   f[0]=(-,-,-) f[1]=(-,-,+) f[2]=(-,+,-) f[3]=(-,+,+)
//   f[4]=(+,-,-) f[5]=(+,-,+) f[6]=(+,+,-) f[7]=(+,+,+)
static real _frac_da_caixa(const real f[8], const real d[3], real fcen) {
    int pos = 0;
    for (int q = 0; q < 8; q++) if (f[q] > 0.0) pos++;
    if (pos == 0) return 0.0;
    if (pos == 8) return 1.0;

    // Parte linear do interpolante trilinear: media e gradiente centrado.
    // O VALOR NO CENTRO, nao a media dos oito cantos.  Num campo linear os
    // dois coincidem, e os oraculos exatos (plano, fatia, tetraedro de canto)
    // nao distinguem.  Num campo CURVO a media dos cantos difere do valor
    // central por (d^2/24) do laplaciano, e esse desvio desloca o plano.
    // Medido no volume da esfera em N=16:  media dos cantos 1.22e-04, valor
    // central 6.11e-05 -- e o valor central bate em quatro digitos com uma
    // segunda implementacao, independente, que usa a distancia exata a'
    // esfera (6.109e-05).  A media dos cantos nao batia com nada.
    real fc = fcen;
    real g[3];
    g[0] = ((f[4]+f[5]+f[6]+f[7]) - (f[0]+f[1]+f[2]+f[3])) / (4.0*d[0]);
    g[1] = ((f[2]+f[3]+f[6]+f[7]) - (f[0]+f[1]+f[4]+f[5])) / (4.0*d[1]);
    g[2] = ((f[1]+f[3]+f[5]+f[7]) - (f[0]+f[2]+f[4]+f[6])) / (4.0*d[2]);

    real gm = sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
    if (gm <= 0.0) return (fc > 0.0) ? 1.0 : 0.0;   // campo plano, sem direcao

    // O LADO DE DENTRO, DIRETO.  Dentro e' func3 > 0, isto e'
    // fc + g.(X - centro) > 0.  Com n = -g/|g| isso vira  n.(X - centro) < s,
    // s = fc/|g|, e o lado de dentro passa a ser o NEGATIVO do plano -- que e'
    // o que a formula fechada calcula.  Importa fazer assim, e nao calcular o
    // lado de fora e subtrair de V: numa celula quase vazia o lado de fora
    // vale quase V, e a subtracao perde todos os digitos da resposta.  Medido:
    // pela subtracao, o tetraedro de canto errava por 100% do seu volume.
    real n[3];
    for (int i = 0; i < 3; i++) n[i] = -g[i]/gm;
    real s = fc/gm;
    real alfa = s;
    for (int i = 0; i < 3; i++) alfa += 0.5*n[i]*d[i];

    // Reflexao dos eixos de normal negativa: x_i -> L_i - x_i leva
    // n_i x_i a n_i L_i + |n_i| x_i', logo alfa ganha |n_i| L_i.  Nao muda o
    // volume, so' o canto de referencia.
    real na[3];
    for (int i = 0; i < 3; i++) {
        na[i] = fabs(n[i]);
        if (n[i] < 0.0) alfa += na[i]*d[i];
    }

    const real V = d[0]*d[1]*d[2];
    real soma = na[0]*d[0] + na[1]*d[1] + na[2]*d[2];
    if (alfa <= 0.0)  return 0.0;
    if (alfa >= soma) return 1.0;

    // Quais eixos sobrevivem.  O limiar: abaixo de 1e-6 a reducao erra O(tol)
    // -- o plano e' reinclinado por esse angulo -- enquanto a formula 3D erra
    // O(eps/prod(n_i)), que ja' passou de 1e-6 ai'.  No caso comum (level-set
    // alinhado a um eixo) as componentes sao ZERO exatos e a reducao e' exata.
    const real tol = 1e-6;
    int viv[3], nv = 0;
    for (int i = 0; i < 3; i++) if (na[i] > tol) viv[nv++] = i;

    real dentro;
    if (nv == 3) {
        dentro = _sz3(na, alfa, d);
    } else if (nv == 2) {
        int a = viv[0], b = viv[1], c = 3 - a - b;   // o eixo que caiu
        dentro = _sz2(na[a], na[b], alfa, d[a], d[b]) * d[c];
    } else if (nv == 1) {
        int a = viv[0];
        dentro = _sz1(na[a], alfa, d[a]) * (V / d[a]);
    } else {
        return (fc > 0.0) ? 1.0 : 0.0;               // normal toda residual
    }
    real fr = dentro / V;
    return (fr < 0.0) ? 0.0 : ((fr > 1.0) ? 1.0 : fr);
}

// Amostra func3 nos oito cantos da caixa de centro `c` e lado `d`.
static void _cantos(const Point c, const real d[3], real f[8]) {
    int q = 0;
    for (int a = -1; a <= 1; a += 2)
        for (int b = -1; b <= 1; b += 2)
            for (int e = -1; e <= 1; e += 2) {
                Point p;
                p[0] = c[0] + 0.5*a*d[0];
                p[1] = c[1] + 0.5*b*d[1];
                p[2] = c[2] + 0.5*e*d[2];
                f[q++] = func3(p);
            }
}

// Fracao de UMA celula, sem subdividir.  Exata se func3 for linear.
real vof3_fracao_celula(Point center, Point delta, real t) {
    real f[8], d[3] = { delta[0], delta[1], delta[2] };
    _cantos(center, d, f);
    return _frac_da_caixa(f, d, func3(center));
}

// Fracao de uma celula com subdivisao NxNxN -- o analogo 3D do laco 16x16
// dos exemplos 2D.  Celula inteiramente de um lado sai na hora; so' a celula
// de interface paga a subdivisao.
real vof3_fracao(Point center, Point delta, real t, int N) {
    real f[8], d[3] = { delta[0], delta[1], delta[2] };
    _cantos(center, d, f);
    int pos = 0;
    for (int q = 0; q < 8; q++) if (f[q] > 0.0) pos++;
    if (pos == 0) return 0.0;
    if (pos == 8) return 1.0;
    if (N < 2) return _frac_da_caixa(f, d, func3(center));

    real ds[3] = { d[0]/N, d[1]/N, d[2]/N };
    real lo[3] = { center[0]-0.5*d[0], center[1]-0.5*d[1], center[2]-0.5*d[2] };
    real acc = 0.0;
    for (int a = 0; a < N; a++)
        for (int b = 0; b < N; b++)
            for (int e = 0; e < N; e++) {
                Point cs;
                cs[0] = lo[0] + (a+0.5)*ds[0];
                cs[1] = lo[1] + (b+0.5)*ds[1];
                cs[2] = lo[2] + (e+0.5)*ds[2];
                real fs[8];
                _cantos(cs, ds, fs);
                acc += _frac_da_caixa(fs, ds, func3(cs));
            }
    return acc / ((real)N*N*N);
}

#endif  // HIGFLOW_EXAMPLES_VOF_GEOMETRY_3D
