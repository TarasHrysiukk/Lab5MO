import java.util.Arrays;

public class Program {
    public static void main(String[] args) {
        MarkwardtMethod mm = new MarkwardtMethod();
        mm.run();
    }
}
class Point {
    // Точка
    private double[] X_N;
    // Градієнт у точці
    private double[] gradientX_N = new double[2];
    // Значення функції у точці
    private double functionX_N;
    // Значення норми градієнту
    private double gradientNorm;
    private static int counter = 0;

    public Point(double[] x_N) {
        this.X_N = x_N;
        functionX_N = Defaults.function(X_N[0], X_N[1]);
        gradientX_N[0] = Defaults.derivativeByX1(x_N[0], x_N[1]);
        gradientX_N[1] = Defaults.derivativeByX2(x_N[0], x_N[1]);
        gradientNorm = Math.sqrt(Math.pow(gradientX_N[0], 2) + Math.pow(gradientX_N[1], 2));
        System.out.println(toString());
        counter++;
    }

    public double[] getX_N() {
        return X_N;
    }

    public double[] getGradientX_N() {
        return gradientX_N;
    }

    public double getFunctionX_N() {
        return functionX_N;
    }

    public double getGradientNorm() {
        return gradientNorm;
    }

    @Override
    public String toString() {
        return "Point "+ counter+ "\t" + Arrays.toString(X_N) +
                "\t\tf(x)=" + functionX_N +
                "\t\tgradient" + Arrays.toString(gradientX_N) +
                "\t\tgradient norm=" + gradientNorm;
    }
}
class Defaults {
    public static final double[] X_0 = {7.0, 7.2};
    public static final double EPSILON = 0.1;
    public static final int MAX_ITERATIONS = 100;
    public static  double LAMBDA = 1;
    public static final double[][] HESSE_MATRIX = {{10,2.8},{2.8,8.4}};//відповідно до коефіційєнтів у градієнті


    public static double function(double x1, double x2) {
        return 5* Math.pow(x1, 2) + 2.8 * x1 * x2 + 4.2 * Math.pow(x2, 2) ;
    }

    public static double derivativeByX1(double x1, double x2) {
        return 10 * x1 + 2.8 * x2 ;
    }

    public static double derivativeByX2(double x1, double x2) {
        return 2.8 * x1 + 8.4 * x2;
    }
}
class MatrixHelper {
    // Додавання матриць
    public double[][] addMatrices(double[][] a, double[][] b) {
        double[][] c = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                c[i][j] = a[i][j] + b[i][j];
            }
        }
        return c;
    }

    public double[] addMatrices(double[] a, double[] b) {
        double[] c = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            c[i] = a[i] + b[i];
        }
        return c;
    }

    // Множення матриці на число
    public double[][] multiplyMatrices(double[][] a, double coefficient) {
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                a[i][j] = a[i][j] * coefficient;
            }
        }
        return a;
    }

    // Множення двох матриць
    public double[][] multiplyMatrices(double A[][], double B[][]) {
        int row1 = A.length;
        int row2 = B.length;
        int col1 = A[0].length;
        int col2 = B[0].length;
        double C[][] = new double[row1][col2];
        if (row2 != col1) {     // Check if multiplication is Possible
            System.out.println("Multiplication Not Possible");
            return null;
        }
        for (int i = 0; i < row1; i++) {
            for (int j = 0; j < col2; j++) {
                for (int k = 0; k < row2; k++)
                    C[i][j] += A[i][k] * B[k][j];
            }
        }
        return C;
    }

    // Переведення n*1 матриці в 1*n
    public double[][] turn1DInto2D(double[] arr) {
        double[][] res = new double[arr.length][1];
        for (int i = 0; i < arr.length; i++) {
            res[i][0] = arr[i];
        }
        return res;
    }

    // Переведення 1*n матриці в n*1
    public double[] turn2DInto1D(double[][] arr) {
        double[] res = new double[arr.length];
        for (int i = 0; i < arr.length; i++) {
            res[i] = arr[i][0];
        }
        return res;
    }

    // Транспонування матриці
    public double[][] transposeMatrix(double[][] m) {
        double[][] temp = new double[m[0].length][m.length];
        for (int i = 0; i < m.length; i++)
            for (int j = 0; j < m[0].length; j++)
                temp[j][i] = m[i][j];
        return temp;
    }

    // Інверсія матриці (піднесення в -1 степінь)
    public double[][] inverseMatrix(double[][] A) {
        double temp;
        int N = A.length;
        double[][] E = new double[N][N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                E[i][j] = 0f;

                if (i == j) {
                    E[i][j] = 1f;
                }
            }
        }
        for (int k = 0; k < N; k++) {
            temp = A[k][k];
            for (int j = 0; j < N; j++) {
                A[k][j] /= temp;
                E[k][j] /= temp;
            }
            for (int i = k + 1; i < N; i++) {
                temp = A[i][k];
                for (int j = 0; j < N; j++) {
                    A[i][j] -= A[k][j] * temp;
                    E[i][j] -= E[k][j] * temp;
                }
            }
        }
        for (int k = N - 1; k > 0; k--) {
            for (int i = k - 1; i >= 0; i--) {
                temp = A[i][k];
                for (int j = 0; j < N; j++) {
                    A[i][j] -= A[k][j] * temp;
                    E[i][j] -= E[k][j] * temp;
                }
            }
        }
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] = E[i][j];
            }
        }
        return E;
    }
}
class MarkwardtMethod {
    private int counter = 0;
    MatrixHelper matrixHelper = new MatrixHelper();

    public void run() {
        // Початкова точка
        Point currentPoint = new Point(Defaults.X_0);
        Point nextPoint;
        while (counter < Defaults.MAX_ITERATIONS) {
            if (currentPoint.getGradientNorm() < Defaults.EPSILON) {
                System.out.println("Gradient norm < EPSILON");
                break;
            }
            nextPoint = new Point(matrixHelper.addMatrices(currentPoint.getX_N(), findS(currentPoint)));
            if (nextPoint.getFunctionX_N() < currentPoint.getFunctionX_N()) {
                Defaults.LAMBDA = Defaults.LAMBDA / 2;
                counter++;
            } else {
                Defaults.LAMBDA = Defaults.LAMBDA * 2;
            }
            currentPoint=nextPoint;
        }
        System.out.println("========================================================");
        System.out.println("Answer is:\t" + currentPoint.toString());
    }

    // Крок вздовж напрямку антиградієнту   S1 = -[H0+λ0I]-1▽f(x0)
    private double[] findS(Point p) {
        double[][] singleMatrix = {{1, 0}, {0, 1}};
        // λ0I
        singleMatrix = matrixHelper.multiplyMatrices(singleMatrix, Defaults.LAMBDA);
        // H0+λ0I
        singleMatrix = matrixHelper.addMatrices(Defaults.HESSE_MATRIX, singleMatrix);
        // -[H0+λ0I]
        singleMatrix = matrixHelper.multiplyMatrices(singleMatrix, -1);
        //-[H0+λ0I]^(-1)
        singleMatrix = matrixHelper.inverseMatrix(singleMatrix);
        // ( -[H0+λ0I]^(-1) )*▽f(x0)
        return matrixHelper.turn2DInto1D(
                matrixHelper.multiplyMatrices(
                        singleMatrix, matrixHelper.turn1DInto2D(p.getGradientX_N())));
    }
}
