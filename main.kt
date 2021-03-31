import java.math.BigDecimal
import java.math.MathContext

fun z1(A: my_Matrix, b: Vector, eps: Double) {
    var counter = 0
    var x = b * BigDecimal(1) // copy b
    while ((x - A * x - b).norm() > BigDecimal(eps)) {
        val new_x = A * x + b
        if (new_x.norm() >= x.norm() + BigDecimal(1))
            counter += 1
        else
            counter = 0
        if (counter > 20) {
            counter = -1
            break
        }
        x = new_x
    }
    if (counter == -1) {
        println(0)
    } else {
        print(x)
    }
}

fun test_z1() {
    // 0.2 0.1 0.3
    // 0.1 0.3 0.2
    // 0.3 0.2 0.1
    var m = my_Matrix(arrayOf(arrayOf(0.2, 0.1, 0.3), arrayOf(0.1, 0.3, 0.2), arrayOf(0.3, 0.2, 0.1)))
    var b = Vector(arrayOf(1.0, 0.0, 1.0))
    z1(m, b, 0.001)
}

fun z2(A: my_Matrix, b: Vector, eps: Double) {
    var counter = 0
    var x = b * BigDecimal(1) // copy b
    var (l, u) = LU(A)
    while ((A * x - b).norm() > BigDecimal(eps)) {
        val new_x = l.obr(b - u * x)
        print(b - u * x)
        print(x)
        print(new_x)
        if (new_x.norm() >= x.norm() + BigDecimal(1))
            counter += 1
        else
            counter = 0
        if (counter > 20) {
            counter = -1
            break
        }

        x = new_x
    }
    if (counter == -1) {
        println(0)
    } else {
        print(x)
    }
}

fun test_z2() {
    // 2.0 0.1 0.3
    // 1.0 3.0 0.2
    // 3.0 2.0 1.0
    var m = my_Matrix(arrayOf(arrayOf(2.0, 0.1, 0.3), arrayOf(1.0, 3.0, 0.2), arrayOf(3.0, 2.0, 1.0)))
    var b = Vector(arrayOf(1.0, 0.0, 1.0))
    z2(m, b, 0.001)
}

fun z3(matrix: Matrix,i: Int,j:Int,phi:Double):Matrix{
    matrix.mult_l(G(0, 0, i, j, phi))
    return matrix
}

fun test_z3() {
    var m11 = my_Matrix(arrayOf(arrayOf(1.0, 2.0), arrayOf(-1.0, 1.0)))
    val phi = 0.33333333
    val i = 0
    val j = 1


    m11.mult_l(G(0, 0, i, j, phi))
    print(m11)
}

fun z4(matrix: Matrix):Pair<Matrix,Matrix> = QR_razl_G(matrix)

fun test_z4() {
    var m4 = my_Matrix(arrayOf(arrayOf(2.0, 0.1, 0.3, 0.5), arrayOf(0.1, 3.0, 0.2,0.2), arrayOf(0.3, 0.2, 1.0,1.0), arrayOf(0.5, 0.2, 1.0, 2.0)))
    val (q,r)=z4(m4.copy())
    println(m4)
    println(q)
    println(r)
    println(q*r)
}


fun z5() {
    val myMatrix = my_Matrix(2, 2)
    val vector = Vector(2)


    myMatrix.mult_l(H(vector))
}

fun z6(matrix: Matrix):Pair<Matrix,Matrix> = QR_razl_H(matrix)

fun test_z6() {
    var m4 = my_Matrix(arrayOf(arrayOf(2.0, 0.1, 0.3, 0.5), arrayOf(0.1, 3.0, 0.2,0.2), arrayOf(0.3, 0.2, 1.0,1.0), arrayOf(0.5, 0.2, 1.0, 2.0)))
    val (q,r)=z6(m4.copy())
    println(m4)
    println(q)
    println(r)
    println(q*r)
}


fun z7(A: my_Matrix, b: complex_vector, eps: Double) {
    val mc = MathContext(A.expeps)
    var counter = 0
    var x = b*b.norm().pow(-1,mc) // copy b
    var lam = complex_namber(BigDecimal(1),BigDecimal(0))
    while ((A * x - x * lam).norm() > BigDecimal(eps)) {
        counter+=1
        val new_x = A * x
        lam = A.norm_with_matrix(x)

        if (counter > 100/eps) {
            counter = -1
            break
        }

        x = new_x*new_x.norm().pow(-1,mc)
    }
    if (counter == -1) {
        println(0)
    } else {
        print(lam)
        print(x)
    }
}

fun test_z7(){
    var m1 = my_Matrix(arrayOf(arrayOf(0.5, 1.0), arrayOf(1.0, 0.5)))
    var cv1 = complex_vector(Vector(arrayOf(1.0, 1.5)),Vector(arrayOf(0.3, -0.5)))  // (real vector), (image vector)
    z7(m1,cv1,0.001)
}

fun z8(matrix: Matrix,eps: Double)=QR_algo(matrix, eps)

fun test_z8(){
    var m4 = my_Matrix(arrayOf(arrayOf(2.0, 0.1, 0.3, 0.5), arrayOf(0.1, 3.0, 0.2,0.2), arrayOf(0.3, 0.2, 1.0,1.0), arrayOf(0.5, 0.2, 1.0, 2.0)))
    val (a,q)=z8(m4.copy(),0.0001)
    println(m4)
    println(a)
    println(q)
    println(q*(a*q.transpose()))
}

fun z9(matrix: Matrix)=free_diag(matrix)

fun test_z9(){
    var m4 = my_Matrix(arrayOf(arrayOf(2.0, 0.1, 0.3, 0.5), arrayOf(0.1, 3.0, 0.2,0.2), arrayOf(0.3, 0.2, 1.0,1.0), arrayOf(0.5, 0.2, 1.0, 2.0)))

    val (a,q)=z9(m4.copy())
    println(m4)
    println(a)
    println(q)
    println(q*(a*(q.transpose())))
}


fun z10(matrix: thee_diag_matrix,eps: Double)=QR_algo_free_diag(matrix,eps)

fun test_z10(){
    var m4 = my_Matrix(arrayOf(arrayOf(2.0, 0.1, 0.3, 0.5), arrayOf(0.1, 3.0, 0.2,0.2), arrayOf(0.3, 0.2, 1.0,1.0), arrayOf(0.5, 0.2, 1.0, 2.0)))

    val (a,q____)=z9(m4.copy())
    val (a1,q)=z10(a,0.0001)
    println(a)
    println(a1)
    println(q)
    println(q*(a1*(q.transpose())))
}
fun z11(matrix: thee_diag_matrix,eps: Double)=QR_algo_free_diag_faster(matrix,eps)

fun test_z11(){
    var m4 = my_Matrix(arrayOf(arrayOf(2.0, 0.1, 0.3, 0.5), arrayOf(0.1, 3.0, 0.2,0.2), arrayOf(0.3, 0.2, 1.0,1.0), arrayOf(0.5, 0.2, 1.0, 2.0)))

    val (a,q____)=z9(m4.copy())
    val (a1,q)=z11(a,0.0001)
    println(a)
    println(a1)
    println(q)
    println(q*(a1*(q.transpose())))
}

fun z12(m1: Matrix,m2: Matrix){
    val A1=super_faste_diag(m1, 0.000001.div(m1.n))
    val A2=super_faste_diag(m2, 0.000001.div(m1.n))
    if ((A1-A2).norm()<BigDecimal(0.000003)){
        print(1)
    }else{
        print(0)
    }
}

fun z13_1(n:Int,d:Int){
    val A1=super_faste_diag(z13_1(n), 0.001.div(d))

    var l =List<BigDecimal>(A1.n) { i -> A1[i, i] }.sortedDescending()

    println(l[1].max(l[A1.n-1]).divide(BigDecimal(d),MathContext(10)))
}
fun z13_2(n:Int,d:Int){
    val A1=super_faste_diag(z13_2(n), 0.001.div(d))

    var l =List<BigDecimal>(A1.n) { i -> A1[i, i] }.sortedDescending()

    println(l[1].max(l[A1.n-1]).divide(BigDecimal(d),MathContext(10)))
}

fun test_z13(){
    z13_1(2,8)
    z13_2(5,3)
}


fun main() {


    var m1 = my_Matrix(arrayOf(arrayOf(0.5, 1.0), arrayOf(1.0, 0.5)))
    var m11 = my_Matrix(arrayOf(arrayOf(1.0, 2.0), arrayOf(-1.0, 1.0)))
    var m4 = my_Matrix(arrayOf(arrayOf(2.0, 0.1, 0.3, 0.5), arrayOf(0.1, 3.0, 0.2,0.2), arrayOf(0.3, 0.2, 1.0,1.0), arrayOf(0.5, 0.2, 1.0, 2.0)))
    var m3 = my_Matrix(arrayOf(arrayOf(2.0, 0.1, 0.3), arrayOf(0.1, 3.0, 0.2), arrayOf(0.3, 0.2, 1.0)))
    var v1 = Vector(arrayOf(1.0, 1.0))
    var cv1 = complex_vector(Vector(arrayOf(1.0, 1.5)),Vector(arrayOf(0.3, -0.5)))
    val u= U(2,2)
    u.copy_From(my_Matrix(arrayOf(arrayOf(2.0, 0.0),arrayOf(2.0, 0.1))))

    // если в тесте много print, и если начальная матрица == конечной, то тест = ок
    test_z1()
    test_z2()
    test_z3()
    test_z4()
    test_z6()
    test_z7()
    test_z8()
    test_z9()
    test_z10()
    test_z11()
    test_z13()

}