// Arcabouco minimo de teste para a HiGTree.
//
// Sem dependencia externa de proposito: a suite ATF de 2020 (higtree/atf-tests)
// parou de rodar porque o `atf-c` nao esta' instalado em lugar nenhum, e 71 casos
// viraram texto morto.  Aqui o unico requisito e' um compilador C.
//
// Cada programa imprime UMA LINHA POR CASO, que o ci/run_higtree_tests.py tabula:
//
//     caso <nome> PASS
//     caso <nome> FAIL <o que falhou, com os numeros>
//
// A linha de falha carrega os valores, nao so' o nome da assercao: um teste que
// diz "esperado 1.234567, obtido 1.230000, tolerancia 1e-09" poupa a rodada de
// depuracao que um "assertion failed" obriga.
#ifndef HIGTREE_TESTING_H
#define HIGTREE_TESTING_H

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int  _t_cases = 0;
static int  _t_failed = 0;
static int  _t_case_failed = 0;
static char _t_case_name[128] = "";

static void _t_case_end(void) {
    if(_t_case_name[0] == '\0') return;
    if(!_t_case_failed) printf("caso %s PASS\n", _t_case_name);
    _t_case_name[0] = '\0';
}

// Abre um caso.  Fecha o anterior, se houver.
static void t_case(const char *nome) {
    _t_case_end();
    snprintf(_t_case_name, sizeof _t_case_name, "%s", nome);
    _t_case_failed = 0;
    _t_cases++;
}

static void _t_fail(const char *fmt, ...) __attribute__((format(printf, 1, 2)));
static void _t_fail(const char *fmt, ...) {
    if(!_t_case_failed) { _t_case_failed = 1; _t_failed++; }
    printf("caso %s FAIL ", _t_case_name);
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    printf("\n");
}

// Encerra o programa.  Devolve o codigo de saida.
static int t_end(void) {
    _t_case_end();
    printf("resumo %d casos, %d falharam\n", _t_cases, _t_failed);
    return _t_failed == 0 ? 0 : 1;
}

#define T_CHECK(cond) \
    do { if(!(cond)) _t_fail("%s:%d  condicao falsa: %s", __FILE__, __LINE__, #cond); } while(0)

#define T_CHECK_MSG(cond, ...) \
    do { if(!(cond)) { _t_fail(__VA_ARGS__); } } while(0)

// Igualdade numerica com tolerancia ABSOLUTA, e a mensagem carrega os tres numeros.
#define T_NEAR(obtido, esperado, tol, oque) \
    do { \
        const double _a = (obtido), _b = (esperado), _t = (tol); \
        if(!(fabs(_a - _b) <= _t)) \
            _t_fail("%s: esperado %.12g, obtido %.12g, erro %.3g > tol %.3g", \
                    (oque), _b, _a, fabs(_a - _b), _t); \
    } while(0)

// Limite SUPERIOR de erro, para quando o valor exato nao e' o criterio e sim a
// ordem de grandeza -- o caso do campo de grau acima do que o esquema reproduz.
#define T_BELOW(valor, limite, oque) \
    do { \
        const double _v = (valor), _l = (limite); \
        if(!(_v <= _l)) \
            _t_fail("%s: %.3g excede o limite %.3g (fator %.1fx)", \
                    (oque), _v, _l, _v / _l); \
    } while(0)

#endif
