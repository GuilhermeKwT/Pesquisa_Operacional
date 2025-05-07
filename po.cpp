#include <iostream>
#include <fstream>
#include <regex>
#include <sstream>
#include <math.h>
#include <map>

using namespace std;

// Faz a leitura do arquivo de entrada e o coloca em uma string com as linhas dividadas por \n
string lerTxt()
{
    ifstream file("input.txt");
    if (!file)
    {
        cout << "num abriu :(" << endl;
    }
    string str, input;
    while (getline(file, str))
    {

        input += str;
        input.push_back('n');
    }
    return input;
}

// Calcula o número de variáveis sem contar as de folga
int calcularNumVariaveis(string input)
{
    int i = 0, num = 0;
    while (input[i] != '=')
    {
        i++;
        cout << i << endl;
        cout << input[i] << endl;
    }
    while (input[i] != 'n')
    {
        cout << i << endl;
        cout << input[i] << endl;
        if (input[i] == 'x')
        {
            num++;
        }
        i++;
    }
    cout << "cabo" << endl << endl;
    return num;
}

/**
 * Calcula o tamanho da matriz A com base no número de variáveis e restrições.
 * @param input     [IN]  A string que possuí o input de onde os coeficientes serão lidos.
 * @param i    [OUT]  O número de variaveis contando as de folg.
 * @param j    [OUT]  O número de restrições.
 */
void calcularIJ(string input, int *i, int *j)
{
    int k = 0;
    (*i) = calcularNumVariaveis(input);
    (*j) = 0;

    while (input[k] != 'n')
    {
        k++;
    }

    while (input[k] != '\0')
    {
        if (input[k] == '>' || input[k] == '<')
        {
            (*i)++;
            (*j)++;
            if (input[k + 1] == '=')
            {
                k++;
            }
        }
        else if (input[k] == '=')
        {
            (*j)++;
        }
        k++;
    }
}

// Faz a alocação dinamica de uma matriz
template <typename T>
T **alocaMatriz(int m, int n)
{
    cout << "c receba" << endl;
    T **M;
    int i;

    M = (T **)malloc(sizeof(T *) * m);
    if (M == NULL)
    {
        printf("Memoria insuficiente.\n");
        exit(1);
    }
    for (i = 0; i < m; i++)
    {
        M[i] = (T *)malloc(sizeof(T) * n);
        if (M[i] == NULL)
        {
            printf("Memoria insuficiente.\n");
            exit(1);
        }
    }
    return M;
}

// Libera a memória de uma matriz dinamicamente alocada
template <typename T>
void liberaMatriz(T **M, int m)
{
    int i;
    for (i = 0; i < m; i++)
        free(M[i]);
    free(M);
}

/**
 * Lê os coeficientes das restrições e as variáveis de folga.
 * @param A     [OUT]  A matriz contendo os coeficientes e as variáveis de folga.
 * @param input     [IN]  A string que possuí o input de onde os coeficientes serão lidos.
 * @param numRestricoes    [IN]  O número de restrições presentes no input.
 * @param numVariaveis    [IN]  O número total de variáveis presentes no input, contando as de folga.
 * @param b    [OUT]  O vetor de recursos.
 * @param c    [OUT]  O vetor de custos.
 */
void lerCoeficientes(double **A, string input, int numRestricoes, int numVariaveis, double *b, double *c)
{
    regex termoRegex("(-?\\d*[.]?\\d*)x(\\d+)");
    regex restricaoRegex("(<=|>=|<|>)");
    regex valRestricao(("\\d+[.]?\\d*"));
    smatch match;
    int linha = 0, k = 0, varFolga = 0;

    string custos;
    while (input[k] != 'n' && input[k] != '\0')
    {
        custos.push_back(input[k]);
        k++;
    }
    k++;

    for(int i = 0; i < numVariaveis; i++){
        c[i] = 0;
    }

    auto it = custos.cbegin();
    while (regex_search(it, custos.cend(), match, termoRegex))
    {
        double coef = (match.str(1).empty() ? 1 : (match.str(1) == "-" ? -1 : stod(match.str(1))));
        int varIndex = stoi(match.str(2)) - 1;
        c[varIndex] = coef;
        it = match.suffix().first;
    }

    // Inicializa com zeros
    for (int i = 0; i < numRestricoes; i++)
    {
        for (int j = 0; j < numVariaveis; j++)
        {
            A[i][j] = 0;
        }
    }

    while (linha < numRestricoes)
    {
        // Separa a restrição
        string restricao;
        while (input[k] != 'n' && input[k] != '\0')
        {
            restricao.push_back(input[k]);
            k++;
        }
        k++; // Pular o '\n'

        auto it = restricao.cbegin();
        while (regex_search(it, restricao.cend(), match, termoRegex))
        {
            double coef = (match.str(1).empty() ? 1 : (match.str(1) == "-" ? -1 : stod(match.str(1))));
            int varIndex = stoi(match.str(2)) - 1;
            A[linha][varIndex] = coef;
            it = match.suffix().first;
        }

        // Separa os valores de b da restrição
        if (regex_search(it, restricao.cend(), match, valRestricao))
        {
            b[linha] = stod(match.str());
        }

        // Detecta o tipo da inequação e define variável de folga
        if (regex_search(restricao, match, restricaoRegex))
        {
            string sinal = match.str(1);
            if (sinal == "<=" || sinal == "<")
            {
                A[linha][calcularNumVariaveis(input) + varFolga] = 1.0;
            }
            else if (sinal == ">=" || sinal == ">")
            {
                A[linha][calcularNumVariaveis(input) + varFolga] = -1.0;
            } // Para '=' não adiciona nada (já está inicializado com zero)
            varFolga++;
        }

        linha++;
    }
}

/**
 * Separa uma submatriz de tamanho n-1 x n-1 da matriz M.
 * @param M     [IN]  Matriz original.
 * @param n     [IN]  Tamanho da matriz original.
 * @param i     [IN]  Linha para ser removida.
 * @param j     [IN]  Coluna para ser removida.
 * @retval **R  Matriz reduzida resultante.
 */
template <typename T>
T **matrizParcial(T **M, int n, int i, int j)
{
    T **R = alocaMatriz<T>(n - 1, n - 1);

    int linha = 0;
    for (int k = 0; k < n; k++)
    {
        if (k == i)
        {
            continue;
        }
        int coluna = 0;
        for (int q = 0; q < n; q++)
        {
            if (q == j)
            {
                continue;
            }
            R[linha][coluna] = M[k][q];
            coluna++;
        }
        linha++;
    }
    return R;
}

/**
 * Realiza o cálculo do determinante de uma matriz quadrada.
 * @param M     [IN]  Matriz que será usada para o cálculo.
 * @param n     [IN]  Tamanho da matriz.
 * @retval double  Resultado do determinante.
 */
double determinante(double **M, int n)
{
    if (n == 1)
    {
        return M[0][0];
    }
    else
    {
        double result = 0, **matrizP;
        int i = 0;
        
        for (int j = 0; j < n; j++)
        {
            matrizP = matrizParcial(M, n, i, j);
            result += pow(-1, i + j) * M[i][j] * determinante(matrizP, n - 1);
            liberaMatriz(matrizP, n - 1);
        }
        
        return result;
    }
}

double **multMatriz(double **A, int lA, int cA, double **B, int lB, int cB){
    double **matrizFinal;
    matrizFinal = alocaMatriz<double>(lA, cB);
    for(int i=0; i<lA; i++){
        for(int j=0; j<cB; j++){
            double somaMult=0;
            for(int k=0; k<cA; k++){
                somaMult = somaMult + (A[i][k] * B[k][j]);
            }
            matrizFinal[i][j] = somaMult;
        }
    }
    return matrizFinal;
}

double **inversa(double **M, int n){
    double **I, **P;
    P = alocaMatriz<double>(n, n);
    I = alocaMatriz<double>(n, n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            P[i][j] = M[i][j]; 
            if (i == j){
                I[i][j] = 1;
            }
            else {
                I[i][j] = 0;
            }
        }
    }
    for (int i = 0; i < n; i++) {
        double divisor = P[i][i];
        for (int j = 0; j <= n; j++){
            P[i][j] /= divisor;
            I[i][j] /= divisor;
        }

        for (int j = 0; j < n; j++){
            if (j != i) {
                double multiplier = P[j][i];
                for (int k = 0; k <= n; k++){
                    P[j][k] -= multiplier * P[i][k];
                    I[j][k] -= multiplier * I[i][k];
                }      
            }
        }
    }
    liberaMatriz(P, n);
    return I;

}

int fatorial(int n){
    if (n == 1){
        return n;
    }
    else {
        return n * fatorial(n - 1);
    }
}

/*
void escolherColunas(double **M, int m, int n, int *B, int *N){
    int conjuntosTestados[fatorial(n)][m], numConjuntos, colunasTemp[m];
    double **MatrizTemp;
    for (int i = 0; i < m; i++){
        colunasTemp[i] = i;
    }
    for (int coluna : colunasTemp){

    }

        B[0] = ((m - 1) % rand()) + 1;
        int aux = B[0];
        bool existe = true;
        for (int i = 1; i < m; i++){
            while (existe){
                aux = ((m - 1) % rand()) + 1;
                int j = 0;
                while (j < i && !existe){
                    if (B[j] == aux){
                        existe = true; 
                    }
                    else {
                        existe = false;
                    }
                    j++;
                }
            }
            B[i] = aux;
        }
    
}
*/

template <typename T>
void imprimeMatriz(T **M, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << '\t' << M[i][j];
        }
        cout << endl;
    }
}

template <typename T>
void impremeVetor(T *v, int n)
{
    for (int i = 0; i < n; i++)
    {
        cout << v[i] << " ";
    }
    cout << endl;
}

int main()
{
    string input = lerTxt();
    int numVariaveis, numRestricoes;
    calcularIJ(input, &numVariaveis, &numRestricoes);

    cout << input << endl;

    double **A = alocaMatriz<double>(numRestricoes, numVariaveis);
    double b[numRestricoes], c[numVariaveis];

    lerCoeficientes(A, input, numRestricoes, numVariaveis, b, c);

    cout << "Vetor c: " << endl;
    impremeVetor(c, numVariaveis);

    cout << endl
         << "Matriz A: " << endl;
    imprimeMatriz(A, numRestricoes, numVariaveis);

    cout << endl
         << "Vetor b: " << endl;
    impremeVetor(b, numRestricoes);

    liberaMatriz(A, numVariaveis);

    /*
    double *porra[3];
    porra = inversa(porra, 3);
    */
    double **M = alocaMatriz<double>(3, 3);
    M[0][0] = 2;
    M[0][1] = 1;
    M[0][2] = 4;
    M[1][0] = 0;
    M[1][1] = 2;
    M[1][2] = 1;
    M[2][0] = 3;
    M[2][1] = 0;
    M[2][2] = 5;
    cout << determinante(M, 3) << endl;
    double **Ml;

    Ml = inversa(M, 3);
    imprimeMatriz(M, 3, 3);

    cout << endl << endl;

    imprimeMatriz(Ml, 3, 3);

    liberaMatriz(Ml, 3);
    cout << "Jean" << endl << endl << endl;
    
    double **I;
    I = alocaMatriz<double>(3, 3);
    cout << "Jean2" << endl << endl << endl;
    I[0][0] = 1;
    I[0][1] = 0;
    I[0][2] = 0;
    I[1][0] = 0;
    I[1][1] = 1;
    I[1][2] = 0;
    I[2][0] = 0;
    I[2][1] = 0;
    I[2][2] = 1;
    
    cout << "biava" << endl;
    M = multMatriz(M, 3, 3, I, 3, 3);
    
    cout << "tashuka" << endl;
    imprimeMatriz(M, 3, 3);

    /*

    double **N = alocaMatriz<double>(2, 2);
    N[0][0] = 1;
    N[0][1] = 5;
    N[1][0] = 2;
    N[1][1] = 3;
    cout << determinante(N, 2) << endl;

    liberaMatriz(A, numRestricoes);
    liberaMatriz(M, 3);
    liberaMatriz(N, 2);
    */
    return 0;
}

