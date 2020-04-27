import sun.invoke.empty.Empty;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

public class Gaussianfields {

    public static void main(String[] args) {
        long start = System.nanoTime();
        parseFile(new File(args[0]), new File(args[1]));
        long end = System.nanoTime();
        System.out.println((end - start) / 1000000);
    }

    private static void parseFile(File inputfile, File outputfile) {

        try {
            Scanner sc = new Scanner(inputfile);

            // This is the exponent of the field
            int F = sc.nextInt();

            // This is the Polynomial that generates the field
            int decimalp = sc.nextInt();

            String binaryp = Integer.toBinaryString(decimalp);

            String[] binaryplist = binaryp.split("");
            int[] p = new int[binaryplist.length];
            for (int i = p.length - 1; i >= 0; i--) {
                p[p.length - 1 - i] = Integer.parseInt(binaryplist[i]);
            }

            //b
            sc.nextLine();
            String aCoeffsAsString = sc.nextLine();
            String[] parameters = aCoeffsAsString.split(", ");
            int m = Integer.parseInt(parameters[0].trim());
            int n = Integer.parseInt(parameters[1].trim());

            int[] b = new int[m];

            // A
            ArrayList<int[]> Atemp = new ArrayList<>();
            int counter = 0;
            while (sc.hasNextLine()) {
                aCoeffsAsString = sc.nextLine().replaceAll("^ +| +$|( )+", "$1");
                String[] coeffsAsList = aCoeffsAsString.split(" ");
                int[] coeffList = new int[coeffsAsList.length - 1];
                for (int i = 0; i < coeffsAsList.length - 1; i++) {
                    coeffList[i] = Integer.parseInt(coeffsAsList[i]);
                }
                Atemp.add(coeffList);
                b[counter] = Integer.parseInt(coeffsAsList[coeffsAsList.length - 1]);
                counter++;
            }
            int[][] A = new int[Atemp.size()][];
            for (int i = 0; i < Atemp.size(); i++) {
                A[i] = Atemp.get(i);
            }

            Polynomial mp = new Polynomial((short) (p.length - 1), p);

            // This is the field root
            short mod = 2;


            Polynomial[] powers = getFieldElements(mp, F, mod);

            int[][] protable;
            int[][] divtable;
            int[] invtable;

            protable = getCalculationTables(powers, mod, F).get(0);
            divtable = getCalculationTables(powers, mod, F).get(1);

            ArrayList<int[]> x = gaussianElimination(A, b, protable, divtable);

            FileWriter writer = new FileWriter(outputfile, false);
            if (x.isEmpty()) {
                for (int i = 0; i < b.length; i++) {
                    writer.write("0");
                    if (i != b.length - 1) {
                        writer.write(" ");
                    }
                }
                writer.close();
            } else if (x.get(0).length > 1) {
                for (int i = 0; i < x.size(); i++) {
                    for (int j = 0; j < x.get(i).length; j++) {
                        writer.write(Integer.toString(x.get(i)[j]));
                        if (j != x.get(i).length - 1) {
                            writer.write(" ");
                        }
                    }
                    if (i != x.size() - 1) {
                        writer.write("\n");
                    }
                }
            } else {
                for (int i = 0; i < x.size(); i++) {
                    writer.write(Integer.toString(x.get(i)[0]));
                    if (i != x.size() - 1) {
                        writer.write(" ");
                    }
                }
            }

            writer.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void combinations(int n, int[] arr, ArrayList<int[]> list) {

        int numArrays = (int) Math.pow(arr.length, n);

        if (numArrays > 10000) {
            numArrays = 10000;
        }
        for (int i = 0; i < numArrays; i++) {
            list.add(new int[n]);
        }
        for (int j = 0; j < n; j++) {
            int period = (int) Math.pow(arr.length, n - j - 1);
            for (int i = 0; i < numArrays; i++) {
                int[] current = list.get(i);
                int index = i / period % arr.length;
                current[j] = arr[index];
            }
        }
    }

    private static ArrayList<int[]> gaussianElimination(int[][] A, int[] b,
                                                        int[][] protable, int[][] divtable) {

        int n = b.length;
        int m = A[0].length;


        ArrayList<ArrayList<Integer>> rows = new ArrayList<>();
        ArrayList<ArrayList<Integer>> columns = new ArrayList<>();

        if (m > n) {
            for (int i = 0; i < m; i++) {
                ArrayList<Integer> a = new ArrayList<>();
                rows.add(a);
                for (int j = 0; j < n; j++) {
                    if (i == 0) {
                        ArrayList<Integer> c = new ArrayList<>();
                        columns.add(c);
                    }

                    if (A[i][j] != 0) {
                        rows.get(i).add(j);
                        columns.get(j).add(i);
                    }
                }
            }
        } else {
            for (int i = 0; i < n; i++) {
                ArrayList<Integer> a = new ArrayList<>();
                rows.add(a);

                for (int j = 0; j < m; j++) {
                    if (i == 0) {
                        ArrayList<Integer> c = new ArrayList<>();
                        columns.add(c);
                    }

                    if (A[i][j] != 0) {
                        columns.get(j).add(i);
                        rows.get(i).add(j);
                    }
                }
            }
        }

        for (int p = 0; p < Math.min(m, n); p++) {

            //checking and switching pivot
            int max = p;
            if (A[p][p] == 0) {

                for (int i = 0; i < columns.get(p).size(); i++) {
                    if (p <= columns.get(p).get(i)) {
                        max = columns.get(p).get(i);
                        break;
                    }
                }
                int[] temp = A[p];
                A[p] = A[max];
                A[max] = temp;

                int t = b[p];
                b[p] = b[max];
                b[max] = t;

                ArrayList<Integer> temporary = rows.get(p);
                rows.set(p, rows.get(max));
                rows.set(max, temporary);

                for (int i = 0; i < columns.get(p).size(); i++) {
                    if (columns.get(p).get(i) == max) {
                        columns.get(p).remove(i);
                        columns.get(p).add(p);
                    }
                }
            }
            // pivot within A and b
            int columncount = columns.get(p).size();
            for (int i = 0; i < columncount; i++) {


                int indexrow = columns.get(p).get(i);
                if (indexrow <= p) {
                    continue;
                }

                int alpha = divtable[A[indexrow][p]][A[p][p]];

                b[indexrow] = b[indexrow] ^ protable[alpha][b[p]];

                int rowcount = rows.get(indexrow).size();
                for (int j = 0; j < rows.get(p).size(); j++) {
                    int indexcol = rows.get(p).get(j);
                    int colcount = columns.get(indexcol).size();
                    A[indexrow][indexcol] = A[indexrow][indexcol] ^ protable[alpha][A[p][indexcol]];
                    if (A[indexrow][indexcol] == 0) {
                        for (int k = 0; k < rowcount; k++) {
                            if (rows.get(indexrow).get(k) == indexcol) {
                                rows.get(indexrow).remove(k);
                                break;
                            }
                        }
                        for (int k = 0; k < colcount; k++) {
                            if (columns.get(indexcol).get(k) == indexrow) {
                                columns.get(indexcol).remove(k);
                                columncount = columns.get(p).size();
                                break;
                            }
                        }
                    } else {
                        if (!rows.get(indexrow).contains(indexcol)) {
                            rows.get(indexrow).add(indexcol);
                        }
                        if (!columns.get(indexcol).contains(indexrow)) {
                            columns.get(indexcol).add(indexrow);
                            columncount = columns.get(p).size();
                        }
                    }
                }
            }
        }

        ArrayList<int[]> Atemp = new ArrayList<>();
        ArrayList<Integer> btemp = new ArrayList<>();

        //removing zero rows
        for (int i = 0; i < n; i++) {
            if (!Arrays.equals(A[i], new int[m])) {
                Atemp.add(A[i]);
                btemp.add(b[i]);
            }
            if (rows.get(i).isEmpty()) {
                rows.remove(i);
            }
        }

        int[][] Aout = new int[Atemp.size()][m];
        int[] bout = new int[btemp.size()];

        for (int i = 0; i < Atemp.size(); i++) {
            Aout[i] = Atemp.get(i);
            bout[i] = btemp.get(i);
        }
        n = bout.length;

        // back substitution

        int[] x = new int[m];

        ArrayList<int[]> xmatrix = new ArrayList<>();

        if (m > n) {
            for (int o = 0; o < protable.length; o++) {
                x = new int[m];
                x[x.length - 1] = o;
                for (int i = n - 1; i >= 0; i--) {
                    int sum = 0;
                    for (int j = 0; j < rows.get(i).size(); j++) {
                        if (rows.get(i).get(j) > i) {
                            sum = sum ^ protable[Aout[i][rows.get(i).get(j)]][x[rows.get(i).get(j)]];
                        }
                    }
                    if (sum == 0) {
                        x[i] = divtable[bout[i]][Aout[i][i]];
                    } else {
                        if (sum == bout[i]) {
                            x[i] = 0;
                        } else {
                            x[i] = divtable[bout[i] ^ sum][Aout[i][i]];
                        }
                    }
                }
                xmatrix.add(x);
            }
            return checkValues(xmatrix, Aout, bout, protable);
        } else {
            boolean check = false;
            for (int i = Math.min(n - 1, m - 1); i >= 0; i--) {
                int sum = 0;
                for (int j = 0; j < rows.get(i).size(); j++) {
                    if (check) {
                        if (rows.get(i).get(j) > i) {
                            sum = sum ^ protable[Aout[i][rows.get(i).get(j)]][x[rows.get(i).get(j)]];
                        }
                    }
                }
                check = true;
                if (sum == 0) {
                    x[i] = divtable[bout[i]][Aout[i][i]];
                } else {
                    if (sum == bout[i]) {
                        x[i] = 0;
                    } else {
                        x[i] = divtable[bout[i] ^ sum][Aout[i][i]];
                    }
                }
            }
            xmatrix.add(x);
            return xmatrix;
        }
    }

    private static ArrayList<int[][]> getCalculationTables(Polynomial[] powers, short mod, int F) {
        int[][] protable = new int[powers.length][powers.length];
        int[][] divtable = new int[powers.length][powers.length];

        for (int i = 0; i < powers.length; i++) {
            for (int j = 0; j < powers.length; j++) {


                int product = multiplyFieldElements(powers, i, j, mod, F);

                protable[i][j] = product;
                divtable[product][i] = j;
            }
        }
        ArrayList<int[][]> tables = new ArrayList<>();
        tables.add(protable);
        tables.add(divtable);

        return tables;
    }

    private static ArrayList<int[]> checkValues(ArrayList<int[]> xmatrix, int[][] A, int[] b, int[][] protable) {
        ArrayList<int[]> solution = new ArrayList<>();

        for (int i = 0; i < xmatrix.size(); i++) {
            boolean checker = true;
            for (int j = 0; j < A.length; j++) {
                int sum = 0;
                for (int k = 0; k < A[0].length; k++) {
                    sum = sum ^ protable[A[j][k]][xmatrix.get(i)[k]];
                }
                if (sum != b[j]) {
                    checker = false;
                }
            }
            if (checker) {
                solution.add(xmatrix.get(i));
            }
        }
        return solution;
    }

    private static Polynomial getSimplePoly(int rank) {
        int[] coeffs = new int[rank + 1];

        coeffs[rank] = 1;
        return new Polynomial((short) rank, coeffs);
    }

    private static Polynomial[] getFieldElements(Polynomial mp, int F, short mod) {
        int elemcount = (int) Math.pow(mod, F);
        int[] pres = new int[mp.coefficients.length - 1];
        System.arraycopy(mp.coefficients, 0, pres, 0, mp.coefficients.length - 1);
        Polynomial respol = new Polynomial(new Integer(mp.degree - 1).shortValue(), pres).removeLeadingZeroes();

        Polynomial alfa0 = new Polynomial((short) -1, new int[]{0});
        Polynomial alfa01 = new Polynomial((short) 0, new int[]{1});

        ArrayList<Polynomial> start = new ArrayList<>();
        start.add(alfa0);
        start.add(alfa01);

        for (int i = 2; i <= F; i++) {
            int[] pol = new int[i];
            pol[i - 1] = 1;
            start.add(new Polynomial((short) (i - 1), pol));
        }

        int counter = 0;
        Polynomial[] powers = new Polynomial[elemcount];

        for (int i = 0; i < mp.degree + 1; i++) {
            powers[i] = start.get(i);
            counter++;
        }


        powers[mp.degree + 1] = respol;
        counter++;

        for (int i = mp.degree + 2; i < elemcount; i++) {
            powers[i] = alfa0;
        }

        for (int i = counter; i < elemcount; i++) {
            powers[i] = powers[i - 1].galoisMultiply(powers[2], mod);
        }

        for (int i = 6; i < powers.length; i++) {
            int rank = powers[i].degree;
            while (rank > F - 1) {
                Polynomial exchangeable = getSimplePoly(rank);
                powers[i] = powers[i].galoisAdd(exchangeable, mod);
                powers[i] = powers[i].galoisAdd(powers[rank + 1], mod);
                rank = powers[i].degree;
            }
        }
        return powers;
    }

    private static int getFieldElement(Polynomial[] field, int element) {

        int fieldelement = 0;
        String[] elementstringbinary = Integer.toBinaryString(element).split("");
        int[] elementintbinary = new int[elementstringbinary.length];

        for (int i = 0; i < elementstringbinary.length; i++) {
            elementintbinary[elementstringbinary.length - 1 - i] = Integer.parseInt(elementstringbinary[i]);
        }

        for (int i = 0; i < field.length; i++) {
            if (Arrays.equals(field[i].coefficients, elementintbinary)) {
                fieldelement = i;
                break;
            }
        }

        return fieldelement;
    }

    private static int getDecimalElement(Polynomial p) {
        StringBuilder sb = new StringBuilder();

        for (int i = p.coefficients.length - 1; i >= 0; i--) {
            sb.append(p.coefficients[i]);
        }

        int decimal = Integer.parseInt(sb.toString(), 2);
        return decimal;
    }

    private static int multiplyFieldElements(Polynomial[] field, int firstelement, int secondelement, short mod, int F) {

        if (firstelement == 0 || secondelement == 0) {
            return 0;
        }
        if (firstelement == 1) {
            return secondelement;
        }
        if (secondelement == 1) {
            return firstelement;
        }

        firstelement = getFieldElement(field, firstelement);
        secondelement = getFieldElement(field, secondelement);

        Polynomial pro = field[firstelement].galoisMultiply(field[secondelement], mod);
        int rank = pro.degree;
        while (rank > F - 1) {
            Polynomial exchangeable = getSimplePoly(rank);
            pro = pro.galoisAdd(exchangeable, mod);
            pro = pro.galoisAdd(field[rank + 1], mod);
            rank = pro.degree;
        }
        return getDecimalElement(pro);
    }

    // From this point on is Polynomial Arithmetic over regular fields
    // This part works, proved in the Coding Theory Course

    private static class Polynomial {

        private short degree;
        private int[] coefficients;

        Polynomial(short degree, int[] coefficients) {
            this.degree = degree;
            this.coefficients = coefficients;
        }

        private Polynomial() {
            this.degree = intAsShort(-1);
            this.coefficients = asIntArray(0);
        }

        short getDegree() {
            return degree;
        }

        int[] getCoefficients() {
            return coefficients;
        }

        private Polynomial galoisAdd(Polynomial b, short modulus) {
            if (this.degree == -1) {
                return b;
            }
            if (b.degree == -1) {
                return this;
            }
            boolean biggerPolynomial = this.degree > b.degree;
            short sumDegree = biggerPolynomial ? this.degree : b.degree;
            short degreeDiff = (intAsShort(Math.abs(this.degree - b.degree)));
            int[] sumCoeffs = new int[sumDegree + 1];
            for (int i = 0; i <= sumDegree - degreeDiff; i++) {
                sumCoeffs[i] = (this.coefficients[i] + b.coefficients[i]) % modulus;
            }
            for (int i = sumDegree - degreeDiff + 1; i <= sumDegree; i++) {
                sumCoeffs[i] = biggerPolynomial ? this.coefficients[i] : b.coefficients[i];
            }
            Polynomial sum = new Polynomial(sumDegree, sumCoeffs);
            return sum.removeLeadingZeroes();
        }

        private Polynomial removeLeadingZeroes() {
            short leadingZeroCount = 0;
            for (int i = this.degree; i >= 0; i--) {
                if (this.coefficients[i] == 0) {
                    leadingZeroCount++;
                } else {
                    break;
                }
            }
            this.degree = (intAsShort(this.getDegree() - leadingZeroCount));
            int[] newCoeffs = new int[this.degree + 1];
            System.arraycopy(this.getCoefficients(), 0, newCoeffs, 0, this.degree + 1);
            this.coefficients = newCoeffs;
            return this;
        }

        Polynomial galoisMultiply(Polynomial b, short mod) {
            if (b.degree == -1 || this.degree == -1) {
                return new Polynomial();
            }
            int sizediff = this.coefficients.length + b.coefficients.length - 1;
            short resultDegree = (intAsShort(this.degree + b.degree));
            int[] prod = new int[sizediff];

            for (int i = 0; i < sizediff; i++) {
                prod[i] = 0;
            }

            for (int i = 0; i < this.coefficients.length; i++) {
                for (int j = 0; j < b.coefficients.length; j++) {
                    prod[i + j] += this.coefficients[i] * b.coefficients[j];
                }
            }
            for (int i = 0; i < prod.length; i++) {
                int result = prod[i] % mod;
                if (result < 0) {
                    result += mod;
                }
                prod[i] = result;
            }
            return new Polynomial(resultDegree, prod);
        }

        private int[] asIntArray(int... coeffs) {
            return coeffs;
        }

        private short intAsShort(int n) {
            return ((short) n);
        }

        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append(this.degree).append("\n");
            if (this.coefficients.length == 0) {
                sb.append("0");
                return sb.toString();
            }
            for (int i = 0; i < this.coefficients.length; i++) {
                sb.append(this.coefficients[i]);
                if (i != this.coefficients.length - 1) {
                    sb.append(" ");
                }
            }
            return sb.toString();
        }
    }
}

