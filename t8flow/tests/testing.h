// Arcabouco minimo de teste do t8Flow.  Vindo do `higtree/tests/testing.h`, e
// a razao de ser igual esta' la': a suite ATF de 2020 parou de rodar porque o
// `atf-c` nao esta' instalado em lugar nenhum, e 71 casos viraram texto morto.
// Aqui o unico requisito e' o compilador que o t8code ja' exige.
//
// Cada programa imprime UMA LINHA POR CASO, que o ci/run_t8flow_tests.py tabula:
//
//     caso <nome> PASS
//     caso <nome> FAIL <o que falhou, com os numeros>
//     resumo <n> casos, <k> falharam
//
// A linha de falha carrega os VALORES, nao so' o nome da assercao: "esperado
// 1.234567, obtido 1.230000, tolerancia 1e-09" poupa a rodada de depuracao que
// um "assertion failed" obriga.
//
// A linha de RESUMO nao e' enfeite: sem ela, "nenhuma falha" e "nenhum caso
// rodou" sao indistinguiveis -- que e' exatamente como se apresentam binario
// obsoleto e segfault no meio da execucao.
#ifndef T8FLOW_TESTING_H
#define T8FLOW_TESTING_H

#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <cmath>

static int  _t_cases = 0;
static int  _t_failed = 0;
static int  _t_case_failed = 0;
static char _t_case_name[128] = "";

static void _t_case_end(void) {
    if (_t_case_name[0] == '\0') return;
    if (!_t_case_failed) printf("caso %s PASS\n", _t_case_name);
    _t_case_name[0] = '\0';
}

//! Abre um caso.  Fecha o anterior, se houver.
static void t_case(const char *nome) {
    _t_case_end();
    snprintf(_t_case_name, sizeof _t_case_name, "%s", nome);
    _t_case_failed = 0;
    _t_cases++;
}

static void _t_fail(const char *fmt, ...) __attribute__((format(printf, 1, 2)));
static void _t_fail(const char *fmt, ...) {
    if (!_t_case_failed) { _t_case_failed = 1; _t_failed++; }
    printf("caso %s FAIL ", _t_case_name);
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    printf("\n");
}

//! Encerra o programa e devolve o codigo de saida.
static int t_end(void) {
    _t_case_end();
    printf("resumo %d casos, %d falharam\n", _t_cases, _t_failed);
    return _t_failed == 0 ? 0 : 1;
}

#define T_CHECK(cond) \
    do { if (!(cond)) _t_fail("%s:%d  condicao falsa: %s", __FILE__, __LINE__, #cond); } while (0)

#define T_CHECK_MSG(cond, ...) \
    do { if (!(cond)) { _t_fail(__VA_ARGS__); } } while (0)

//! Igualdade numerica com tolerancia ABSOLUTA; a mensagem carrega os tres numeros.
#define T_NEAR(obtido, esperado, tol, oque) \
    do { \
        const double _a = (obtido), _b = (esperado), _t = (tol); \
        if (!(std::fabs(_a - _b) <= _t)) \
            _t_fail("%s: esperado %.12g, obtido %.12g, erro %.3g > tol %.3g", \
                    (oque), _b, _a, std::fabs(_a - _b), _t); \
    } while (0)

//! Limite SUPERIOR de erro, para quando o criterio e' ordem de grandeza e nao
//! valor exato -- o caso do campo de grau acima do que o esquema reproduz.
#define T_BELOW(valor, limite, oque) \
    do { \
        const double _v = (valor), _l = (limite); \
        if (!(_v <= _l)) \
            _t_fail("%s: %.3g excede o limite %.3g (fator %.1fx)", \
                    (oque), _v, _l, _v / _l); \
    } while (0)

//! CONTRASTE: exige que `pior` seja ao menos `fator` vezes maior que `melhor`.
//! Sem isto, "o erro e' pequeno" nao distingue reproducao exata de acerto por
//! acidente -- e foi assim que a exatidao polinomial escondeu parede errada
//! com campo linear na suite da HiGTree.
#define T_CONTRASTE(melhor, pior, fator, oque) \
    do { \
        const double _m = (melhor), _p = (pior), _f = (fator); \
        if (!(_p > _m * _f)) \
            _t_fail("%s: esperava contraste > %.3gx, obtive %.3gx " \
                    "(melhor %.3g, pior %.3g) -- o oraculo nao discrimina", \
                    (oque), _f, (_m > 0 ? _p / _m : 0.0), _m, _p); \
    } while (0)

#endif
