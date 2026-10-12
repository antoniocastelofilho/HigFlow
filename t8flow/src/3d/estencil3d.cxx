// Implementacao do construtor de estencil.  O PORQUE da forma -- `rhs` separado,
// parada ao fechar o anel, suporte local -- esta' no `estencil3d.hxx`.  Aqui ficam
// as decisoes de implementacao:
//
//   1. os aneis sao construidos por VIZINHANCA DE FACE, nao por raio geometrico.
//      Raio ignora a topologia: numa transicao 2:1 ele pega celulas que nao sao
//      vizinhas e perde as que sao;
//   2. o anel e' SEMPRE fechado por inteiro antes de decidir parar.  Decidir no
//      meio faz o suporte depender da ordem de visita;
//   3. o ponto de CONTORNO entra no suporte com id -1 e com o seu valor de
//      Dirichlet guardado a parte.  O peso dele multiplica o valor e vai para o
//      `rhs`; ele NAO aparece em `ids`;
//   4. o conjunto de ja'-vistos e' por ID, nao por posicao: comparar posicao em
//      ponto flutuante para decidir "ja' tenho esta celula" e' como se duplica
//      celula perto de interface de bloco.
#include "estencil3d.hxx"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>

namespace t8flow3d {

ConstrutorEstencil::ConstrutorEstencil(const Malha3 &m, int ordem, int max_pontos,
                                       int max_aneis)
    : m_(&m), ordem_(ordem), max_pontos_(max_pontos), max_aneis_(max_aneis),
      mls_(ordem, max_pontos)
{
  if (max_aneis < 1)
    throw std::invalid_argument("ConstrutorEstencil: max_aneis tem de ser >= 1");
}

namespace {

//! (bloco, folha) de um id global -- a malha indexa por par, o vizinho vem por id.
//
// Varredura: O(folhas).  E' aceitavel aqui porque a selecao de suporte ja' e' um
// caminho de montagem, nao de laco quente; quando for, isto vira um indice.  Dito
// para nao ser suposto rapido.

}  // namespace

ConstrutorEstencil::Suporte ConstrutorEstencil::suporte_em(const Ponto3 &x) const
{
  return (estrategia_ == Suporte3::Caixa) ? suporte_caixa(x) : suporte_aneis(x);
}

//! SUPORTE POR CAIXA GEOMETRICA: os pontos em `c + (i,j,k)*h`.
//
// POR QUE ELE EXISTE, e a medida e' a razao inteira: com o suporte por aneis de
// face, 504 de 512 ajustes de ordem 3 REBAIXARAM a ordem -- 98,4%.  Os pontos
// bastavam em numero (24,4 em media para 20 monomios) e nao em DISTRIBUICAO: a
// cruz tem quase tudo alinhado aos eixos, e os monomios mistos `xy`, `xz`, `yz`
// e `xyz` ficam mal determinados nela.
//
// A CAIXA entrega 27 pontos num cubo 3x3x3, que e' a configuracao classica para
// ajuste tridimensional, e condiciona a matriz normal muito melhor.
//
// SO' E' VIAVEL PORQUE `id_do_ponto` E' O(niveis).  Com a varredura linear que o
// lado 2D teve ate' ser consertada, 27 buscas por ponto de consulta tornariam a
// montagem inviavel -- e a montagem ja' foi quadratica uma vez neste projeto.
//
// O PASSO E' O DA CELULA DO ALVO, e em malha adaptada isso e' uma escolha, nao
// uma obviedade: um vizinho mais fino cai dentro do mesmo passo e entra uma vez
// so' (os ids repetidos sao filtrados); um mais grosso e' alcancado por varios
// pontos e tambem entra uma vez.  O suporte fica geometricamente regular mesmo
// onde a malha nao e'.
ConstrutorEstencil::Suporte ConstrutorEstencil::suporte_caixa(const Ponto3 &x) const
{
  const long id0 = m_->id_do_ponto(x);
  if (id0 < 0)
    throw std::out_of_range("ConstrutorEstencil: ponto fora do dominio");
  const long x0 = m_->indice_do_id(id0);
  if (x0 < 0)
    throw std::runtime_error("ConstrutorEstencil: id " + std::to_string(id0) +
                             " nao esta' nas folhas locais nem na franja");

  const int minimo = Mls3::min_pontos(ordem_);

  // O RAIO MINIMO NAO SAI DA CONTAGEM DE PONTOS, E SIM DA ALGEBRA.
  //
  // Uma caixa de raio `r` tem `2r+1` valores DISTINTOS de coordenada por
  // direcao, e um polinomio de grau `g` precisa de `g+1` valores distintos para
  // ser determinado naquela direcao.  Com menos, os monomios ficam linearmente
  // DEPENDENTES nos dados: num cubo 3x3x3, onde `x` so' assume `-h`, `0` e `h`,
  // vale `x^3 = h^2 * x` EXATAMENTE -- e a matriz normal e' singular.
  //
  // MEDIDO, e e' o que trouxe esta formula: parando pela contagem, a caixa de
  // raio 1 dava 27 pontos para 20 monomios -- numero de sobra -- e 216 dos 512
  // ajustes de ordem 3 rebaixavam.  Os 216 eram exatamente `6^3`, as celulas do
  // INTERIOR: as de borda cresciam para o raio 2 por falta de pontos e por isso
  // nao rebaixavam.  O sintoma apontou a causa.
  //
  //     ordem 2 -> 3 valores -> raio 1
  //     ordem 3 -> 4 valores -> raio 2
  //     ordem 4 -> 5 valores -> raio 2
  const int raio_minimo = (ordem_ + 1) / 2;

  const Ponto3 c0 = m_->centroide_do_indice(x0);
  const Ponto3 h0 = m_->lados_do_indice(x0);

  Suporte s;
  std::set<long> vistos;

  auto inclui = [&](long id, long idx) {
    if (!vistos.insert(id).second) return;
    s.pontos.push_back(m_->centroide_do_indice(idx));
    s.ids.push_back(id);
    s.a.push_back(1.0);
    s.b.push_back(0.0);
  };
  inclui(id0, x0);

  for (int r = 1; r <= max_aneis_; r++) {
    // A CASCA DE RAIO `r`, e nao a caixa inteira: os de dentro ja' entraram.
    // Percorrer a caixa toda a cada raio repetiria trabalho e, pior, faria o
    // suporte depender da ordem de visita quando o criterio de parada
    // interrompesse no meio.
    for (int i = -r; i <= r; i++)
      for (int j = -r; j <= r; j++)
        for (int k = -r; k <= r; k++) {
          if (std::abs(i) != r && std::abs(j) != r && std::abs(k) != r) continue;
          const Ponto3 p{c0[0] + i * h0[0], c0[1] + j * h0[1], c0[2] + k * h0[2]};
          const long id = m_->id_do_ponto(p);
          if (id < 0) {
            // AS DUAS RESPOSTAS NEGATIVAS SAO DIFERENTES, e confundi-las esconde
            // a segunda -- que e' a que importa.
            //
            // FORA DO DOMINIO e' legitimo: o fechamento de contorno trata, pelo
            // laco de faces abaixo (que sabe POR QUAL face o ponto saiu; a caixa
            // nao sabe, e inventar daria a condicao errada numa quina).
            //
            // DENTRO do dominio e nao achado significa que falta FRANJA, e isso
            // tem de ser CONTADO: um suporte truncado da' pesos diferentes dos
            // do serial, logo um resultado que depende da particao.
            //
            // Sem esta distincao a medida de "quantas camadas bastam" dava ZERO
            // truncamento com qualquer numero de camadas -- um resultado bom
            // demais, que e' sintoma e nao conquista.
            if (m_->dentro_do_dominio(p)) s.truncado_pela_franja = true;
            continue;
          }
          const long idx = m_->indice_do_id(id);
          if (idx < 0) continue;
          inclui(id, idx);
        }

    // A CASCA E' FECHADA POR INTEIRO antes de decidir parar -- parar no meio
    // faria o suporte, e portanto os pesos, depender da ordem de visita.
    //
    // E AS DUAS CONDICOES TEM DE VALER: pontos suficientes E raio suficiente.
    // A contagem sozinha deixa passar o suporte com monomios dependentes; o
    // raio sozinho nao garante numero em malha adaptada, onde celulas grossas
    // vizinhas fazem varios pontos caírem na mesma folha.
    if (r >= raio_minimo && static_cast<int>(s.pontos.size()) >= minimo) break;
  }

  // O FECHAMENTO DE CONTORNO vem pelas FACES da celula do alvo, que e' onde a
  // condicao esta' definida.  A caixa nao sabe por qual face um ponto saiu do
  // dominio; a face sabe.
  if (!regioes_.vazio()) {
    for (int fa = 0; fa < 6; fa++) {
      if (!m_->vizinhos_do_indice(x0, fa).empty()) continue;
      if (!m_->expansivel_do_indice(x0)) continue;   // borda da franja, nao do dominio

      const int d = fa / 2;
      const double sinal = (fa % 2 == 0) ? -1.0 : 1.0;
      Ponto3 parede = c0;
      parede[d] += sinal * h0[d] / 2.0;

      const Contorno &cond = regioes_.em(parede);
      if (cond.tipo == TipoContorno::Dirichlet) {
        s.pontos.push_back(parede);
        s.ids.push_back(-1);
        s.a.push_back(0.0);
        s.b.push_back(cond.g(parede));
        s.paredes++;
      } else {
        Ponto3 fantasma = c0;
        fantasma[d] += sinal * h0[d];
        double A = 0.0, B = 0.0;
        cond.afim(parede, h0[d], &A, &B);
        s.pontos.push_back(fantasma);
        s.ids.push_back(id0);
        s.a.push_back(A);
        s.b.push_back(B);
        s.fantasmas++;
      }
    }
  }

  if (static_cast<int>(s.pontos.size()) < minimo)
    throw std::runtime_error(
        "ConstrutorEstencil: suporte de " + std::to_string(s.pontos.size()) +
        " pontos em " + std::to_string(max_aneis_) + " cascas, minimo " +
        std::to_string(minimo) + " para ordem " + std::to_string(ordem_));
  if (s.truncado_pela_franja) truncados_++;
  return s;
}

ConstrutorEstencil::Suporte ConstrutorEstencil::suporte_aneis(const Ponto3 &x) const
{
  const long id0 = m_->id_do_ponto(x);
  if (id0 < 0)
    throw std::out_of_range("ConstrutorEstencil: ponto fora do dominio");

  const int minimo = Mls3::min_pontos(ordem_);

  Suporte s;
  std::set<long> vistos;        // por ID, nao por posicao -- ver o cabecalho
  std::vector<long> fronteira;  // os ids do anel corrente

  const long x0 = m_->indice_do_id(id0);
  if (x0 < 0)
    throw std::runtime_error("ConstrutorEstencil: id " + std::to_string(id0) +
                             " nao esta' nas folhas locais nem na franja");

  vistos.insert(id0);
  fronteira.push_back(id0);
  s.pontos.push_back(m_->centroide_do_indice(x0));
  s.ids.push_back(id0);
  s.a.push_back(1.0);
  s.b.push_back(0.0);

  for (int anel = 1; anel <= max_aneis_; anel++) {
    // O ANEL E' FECHADO POR INTEIRO antes de qualquer decisao de parada.  Parar no
    // meio faria o suporte -- e portanto os pesos -- depender da ordem de visita.
    std::vector<long> proximo;
    for (long id : fronteira) {
      const long xi = m_->indice_do_id(id);
      if (xi < 0) { s.truncado_pela_franja = true; continue; }

      // A FOLHA DA ULTIMA CAMADA DA FRANJA NAO SE EXPANDE, e distinguir isso de
      // "face no contorno do dominio" e' obrigatorio: as duas situacoes dao lista
      // de vizinhos VAZIA, e confundi-las faria o esquema aplicar condicao de
      // contorno no MEIO do dominio -- em paralelo, ao longo de toda a fronteira
      // de particao, e em silencio.
      //
      // Aqui o suporte pode ficar TRUNCADO, e o contador existe para que isso
      // nao passe calado: um suporte truncado da' pesos diferentes dos do
      // serial, logo um resultado que DEPENDE do numero de processos.  Com
      // `exige_camadas(aneis)` isso nao deve acontecer -- e o contador e' como se
      // afirma que nao aconteceu.
      if (!m_->expansivel_do_indice(xi)) { s.truncado_pela_franja = true; continue; }

      const Ponto3 c_xi = m_->centroide_do_indice(xi);
      const Ponto3 h_xi = m_->lados_do_indice(xi);
      for (int fa = 0; fa < 6; fa++) {
        const std::vector<long> vz = m_->vizinhos_do_indice(xi, fa);
        if (vz.empty()) {
          // FACE DE CONTORNO.  Sem condicao instalada o contorno nao contribui, e
          // o suporte vai faltar pontos -- que e' melhor que contribuir com zero.
          if (regioes_.vazio()) continue;

          const Ponto3 c = c_xi;
          const Ponto3 h = h_xi;
          const int d = fa / 2;
          const double sinal = (fa % 2 == 0) ? -1.0 : 1.0;

          Ponto3 parede = c;
          parede[d] += sinal * h[d] / 2.0;          // centro da faceta

          // A REGIAO decide qual condicao vale NESTA faceta.  Ponto descoberto ou
          // em duas regioes faz `em()` lancar -- ver o `contorno3d.hxx`.
          const Contorno &cond = regioes_.em(parede);

          if (cond.tipo == TipoContorno::Dirichlet) {
            // DIRICHLET amostra a PAREDE: o valor ali e' conhecido e nao depende de
            // incognita nenhuma.  E' o caminho mais direto, e nao se troca por
            // fantasma de proposito -- refletir (A=-1) amplifica o erro do interior,
            // amostrar nao.
            s.pontos.push_back(parede);
            s.ids.push_back(-1);
            s.a.push_back(0.0);
            s.b.push_back(cond.g(parede));
            s.paredes++;
          } else {
            // NEUMANN e ROBIN usam FANTASMA: um ponto do lado de FORA, espelhado
            // pela faceta, cujo valor e' afim na incognita desta celula.  A
            // distancia entre os centros e' o lado da celula -- dois meios-passos.
            Ponto3 fantasma = c;
            fantasma[d] += sinal * h[d];
            double A = 0.0, B = 0.0;
            cond.afim(parede, h[d], &A, &B);
            s.pontos.push_back(fantasma);
            s.ids.push_back(id);                    // a celula ESPELHADA
            s.a.push_back(A);
            s.b.push_back(B);
            s.fantasmas++;
          }
          continue;
        }
        for (long v : vz) {
          if (vistos.count(v) != 0) continue;
          // O VIZINHO PODE ESTAR NA FRANJA, e entra no suporte do mesmo jeito:
          // o que a franja da' -- centro e id global -- e' tudo que o ajuste
          // precisa de um ponto de suporte.  A versao anterior usava uma busca
          // que so' cobria as folhas locais e seguia com `continue`, o que em
          // paralelo DESCARTARIA o vizinho sem dizer nada.
          const long vx = m_->indice_do_id(v);
          if (vx < 0) { s.truncado_pela_franja = true; continue; }
          vistos.insert(v);
          proximo.push_back(v);
          s.pontos.push_back(m_->centroide_do_indice(vx));
          s.ids.push_back(v);
          s.a.push_back(1.0);
          s.b.push_back(0.0);
        }
      }
    }
    s.aneis = anel;
    if (proximo.empty()) break;            // a malha acabou (dominio pequeno)
    fronteira = proximo;
    if (static_cast<int>(s.pontos.size()) >= minimo) break;   // anel FECHADO, agora para
  }

  // TETO: o ajuste nao aceita mais que `max_pontos`, e cortar o excesso tem de ser
  // por DISTANCIA ao alvo -- cortar pela ordem de insercao devolveria um suporte
  // que depende do percurso.
  if (static_cast<int>(s.pontos.size()) > max_pontos_) {
    std::vector<int> ordem(s.pontos.size());
    for (size_t k = 0; k < ordem.size(); k++) ordem[k] = static_cast<int>(k);
    const auto d2 = [&](int k) {
      const double dx = s.pontos[k][0] - x[0], dy = s.pontos[k][1] - x[1];
      return dx * dx + dy * dy;
    };
    std::stable_sort(ordem.begin(), ordem.end(),
                     [&](int a, int b) { return d2(a) < d2(b); });
    Suporte t;
    t.aneis = s.aneis;
    for (int k = 0; k < max_pontos_; k++) {
      t.pontos.push_back(s.pontos[ordem[k]]);
      t.ids.push_back(s.ids[ordem[k]]);
      t.a.push_back(s.a[ordem[k]]);
      t.b.push_back(s.b[ordem[k]]);
      if (s.ids[ordem[k]] < 0) t.paredes++;
      else if (s.a[ordem[k]] != 1.0 || s.b[ordem[k]] != 0.0) t.fantasmas++;
    }
    s = t;
  }

  if (static_cast<int>(s.pontos.size()) < minimo)
    throw std::runtime_error(
        "ConstrutorEstencil: suporte de " + std::to_string(s.pontos.size()) +
        " pontos em " + std::to_string(max_aneis_) + " aneis, minimo " +
        std::to_string(minimo) + " para ordem " + std::to_string(ordem_) +
        (regioes_.vazio() ? " (nenhuma condicao de contorno instalada: o contorno"
                            " nao contribuiu)" : ""));

  for (size_t a = 0; a < s.pontos.size(); a++)
    for (size_t b = a + 1; b < s.pontos.size(); b++) {
      const double dx = s.pontos[a][0] - s.pontos[b][0];
      const double dy = s.pontos[a][1] - s.pontos[b][1];
      s.diametro = std::fmax(s.diametro, std::sqrt(dx * dx + dy * dy));
    }
  if (s.truncado_pela_franja) truncados_++;
  return s;
}

Estencil ConstrutorEstencil::em(const Ponto3 &x) const
{
  const Suporte s = suporte_em(x);
  const std::vector<double> w = mls_.pesos(x, s.pontos);

  // O `rhs` RECEBE a contribuicao dos pontos de contorno, e eles NAO entram em
  // `ids`.  Quem consome o estencil monta a linha do sistema com `ids`/`pesos` e
  // soma `rhs` ao lado direito -- sem precisar saber qual id era parede.
  return monta(s, w);
}

Estencil ConstrutorEstencil::monta(const Suporte &s, const std::vector<double> &w)
{
  // UM LACO PARA OS TRES TIPOS: `valor_no_ponto = a*u[id] + b`, logo o peso `w`
  // contribui `a*w` no id e `b*w` no `rhs`.  Celula comum tem (a,b)=(1,0), parede
  // de Dirichlet tem id=-1 e (0,g), fantasma tem o id da celula ESPELHADA e os
  // (A,B) da condicao -- e e' por isso que o mesmo id pode aparecer duas vezes:
  // a celula e o fantasma dela.  Os pesos se SOMAM nesse caso, e somar e' o
  // correto; sobrescrever perderia a contribuicao do fantasma em silencio.
  Estencil e;
  e.rhs = 0.0;
  for (size_t k = 0; k < w.size(); k++) {
    e.rhs += s.b[k] * w[k];
    if (s.ids[k] < 0) continue;
    const double peso = s.a[k] * w[k];
    bool somado = false;
    for (size_t j = 0; j < e.ids.size(); j++)
      if (e.ids[j] == s.ids[k]) { e.pesos[j] += peso; somado = true; break; }
    if (!somado) { e.ids.push_back(s.ids[k]); e.pesos.push_back(peso); }
  }
  return e;
}

Estencil ConstrutorEstencil::derivada_em(const Ponto3 &x, int dir) const
{
  const Suporte s = suporte_em(x);
  // A MONTAGEM E' A MESMA do valor -- `(a, b, peso)` -> `(ids, pesos, rhs)` --, e
  // reusa-la aqui nao e' economia: e' o que garante que o fechamento de contorno
  // entre na derivada do mesmo jeito que entra no valor.  Duas montagens
  // divergiriam no `rhs`, e divergir ali so' aparece perto da parede.
  return monta(s, mls_.pesos_derivada(x, s.pontos, dir));
}

Estencil ConstrutorEstencil::laplaciano_em(const Ponto3 &x) const
{
  const Suporte s = suporte_em(x);
  return monta(s, mls_.pesos_laplaciano(x, s.pontos));
}


}  // namespace t8flow3d
