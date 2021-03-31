import java.math.BigDecimal
import java.math.MathContext
import kotlin.math.min
import kotlin.math.sqrt
import kotlin.test.assertEquals


fun LU(matrix: my_Matrix): Pair<L, my_Matrix> {
    val l = L(matrix.n, matrix.m)
    val u = my_Matrix(matrix.n, matrix.m)
    for (j in 0 until matrix.m) {
        for (i in j until matrix.n) {
            l.set(i, j, matrix.get(i, j))
        }
    }
    for (j in 0 until matrix.m) {
        for (i in 0 until j) {
            u.set(i, j, matrix.get(i, j))
        }
    }
    return Pair(l, u)
}

fun QR_razl_G(myMatrix: Matrix): Pair<Matrix, Matrix> {
    var A = myMatrix.copy()
    var Q = E(A.n, A.m)
    val mc = MathContext(A.expeps)
    for (k in 0 until A.m) {
        var i = k
        var v_i = A[i, k]
        while ((v_i < BigDecimal(2E-10)) && (i < A.n)) {
            i++
        }
        if (i == A.n) continue
        for (j in (i + 1) until A.n) {
            v_i = A[i, k]
            val v_j = A[j, k]
            if (v_i.abs() < BigDecimal(2E-10)) continue

            val sum_2 = (v_i.pow(2, mc) + v_j.pow(2, mc)).sqrt(mc)
            if (sum_2.abs() < BigDecimal(2E-10)) continue
//            println(Pair(v_i,sum_2))
            val c = v_i.divide(sum_2, mc)
            val s = -v_j.divide(sum_2, mc)
            var g = G(0, 0, j, i, c, s)

            A.mult_l(g)
            Q.mult_l(g)
        }
        if (i == k) continue

        val pi = 3.141592653589
        var g = G(0, 0, k, i, pi / 2)

        A.mult_l(g)
        Q.mult_l(g)
    }
    return Pair(Q.transpose(), A)
}

fun QR_razl_H(myMatrix: Matrix): Pair<Matrix, Matrix> {
    var A = myMatrix.copy()
    var Q = E(A.n, A.m)

    val mc = MathContext(A.expeps)
    for (i in 0 until A.m) {
        var A2 = sub_matrix(zdvig_n = i, zdvig_m = i, v_matrix = A)
        var Q2 = sub_matrix(zdvig_n = i, zdvig_m = i, v_matrix = Q)
        val e_i = Vector(A2.n)
        e_i[0] = BigDecimal(1)

        val v = A2.gett(0)
        if (v.norm() == BigDecimal(0)) continue
        val u = v * (v.norm().pow(-1, mc))
        if (u == e_i) continue
        val u2 = (u - e_i) * ((u - e_i).norm().pow(-1, mc))

        val h = H(u2)

        A2.mult_l(h)
        Q2.mult_l(h)

    }
    return Pair(Q.transpose(), A)
}

fun QR_algo(matrix: Matrix, eps: Double): Pair<Matrix, Matrix> {
    var (Q, R) = QR_razl_G(matrix)
    var Q_out = Q.copy()
    var A = R * Q
    while (A.norm_ger() > BigDecimal(eps)) {
        val (Q, R) = QR_razl_G(A)
        A = R * Q
        Q_out = Q_out * Q
    }
    return Pair(A, Q_out)
}

fun free_diag(A: Matrix): Pair<thee_diag_matrix, Matrix> {
    val (AA, QQ) = free_diag(A, 0)
    val A_out = thee_diag_matrix(A.n)
    A_out.copy_From(AA)
    return Pair(A_out, QQ)
}

fun free_diag(A: Matrix, zdvig: Int): Pair<Matrix, Matrix> {
    if (A.n - zdvig <= 2) return Pair(A, E(A.n, A.n))

    val mc = MathContext(A.expeps)
    var A1 = sub_matrix(zdvig_n = zdvig, zdvig_m = zdvig, v_matrix = A)
    val e_1 = Vector(A1.n)
    e_1[1] = BigDecimal(1)
    var u1 = A1.gett(0).copy()
    u1[0] = BigDecimal(0)
    u1 = u1 * u1.norm().pow(-1, mc)
    var u = (u1 - e_1) * (u1 - e_1).norm().pow(-1, mc)
    var h = H2(u)
    A1.mult_lr(h)


    var (AA, QQ) = free_diag(A, zdvig + 1)
    var QQ1 = sub_matrix(zdvig_n = zdvig, zdvig_m = zdvig, v_matrix = QQ)
//    QQ1.copy_From(H(u) * QQ1)
//    QQ1.copy_From()
    QQ1.mult_l(H(u))

    return Pair(AA, QQ)
}


fun QR_razl_G_free_diag(mythee_diag_matrix: thee_diag_matrix): Pair<Matrix, Matrix> {
    var A = my_Matrix(mythee_diag_matrix.n, mythee_diag_matrix.m)
    A.copy_From(mythee_diag_matrix)
    var Q = E(A.n, A.m) //E
    val mc = MathContext(A.expeps)
    for (k in 0 until A.m) {
        var i = k
        val v_i = A[i, k]
        while ((v_i == BigDecimal(0)) && (i < A.n) && (i < k + 3)) {
            i++
        }
        if ((i == A.n) || (i == k + 3)) continue
        for (j in (i + 1) until min(A.n, k + 2)) {
            val v_j = A[j, k]
            if (v_i == BigDecimal(0)) continue

            val sum_2 = (v_i.pow(2, mc) + v_j.pow(2, mc)).sqrt(mc)
            val c = v_i.divide(sum_2, mc)
            val s = -v_j.divide(sum_2, mc)
            var g = G(0, 0, j, i, c, s)

            A.mult_l(g)
            Q.mult_l(g)
        }
        if (i == k) continue

        val pi = 3.141592653589
        var g = G(0, 0, k, i, pi / 2)

        A.mult_l(g)
        Q.mult_l(g)
    }

    return Pair(Q.transpose(), A)
}
fun QR_razl_G_free_diag2(mythee_diag_matrix: thee_diag_matrix,Q_out : Matrix): Pair<Matrix, Matrix> {
    var A = my_Matrix(mythee_diag_matrix.n, mythee_diag_matrix.m)
    A.copy_From(mythee_diag_matrix)
    var Q = E(A.n, A.m) //E
    val mc = MathContext(A.expeps)
    for (k in 0 until A.m) {
        var i = k
        val v_i = A[i, k]
        while ((v_i == BigDecimal(0)) && (i < A.n) && (i < k + 3)) {
            i++
        }
        if ((i == A.n) || (i == k + 3)) continue
        for (j in (i + 1) until min(A.n, k + 2)) {
            val v_j = A[j, k]
            if (v_i == BigDecimal(0)) continue

            val sum_2 = (v_i.pow(2, mc) + v_j.pow(2, mc)).sqrt(mc)
            val c = v_i.divide(sum_2, mc)
            val s = -v_j.divide(sum_2, mc)
            var g = G(0, 0, j, i, c, s)

            A.mult_l(g)
            Q.mult_l(g)
            Q_out.mult_l(g)
        }
        if (i == k) continue

        val pi = 3.141592653589
        var g = G(0, 0, k, i, pi / 2)

        A.mult_l(g)
        Q.mult_l(g)
    }

    return Pair(Q.transpose(), A)
}

fun QR_algo_free_diag(matrix: thee_diag_matrix, eps: Double): Pair<thee_diag_matrix, Matrix> {
    var (Q, R) = QR_razl_G_free_diag(matrix)
//    var Q_out = Q.copy()
    var Q_out = my_Matrix(Q.n, Q.m)
    Q_out.copy_From(Q.transpose())
    var A = thee_diag_matrix(R.n)
    A.copy_From(R.mult_r_for_3_diag(Q))
    while (A.norm_ger() > BigDecimal(eps)) {
        val (Q_, R_) = QR_razl_G_free_diag2(A, Q_out)
        A.copy_From(R_.mult_r_for_3_diag(Q_))
//        Q_out = Q_out * Q
//        Q_out.copy_From(Q_out * Q_)

    }

    return Pair(A, Q_out.transpose())
}

fun take_s(mat: Matrix): BigDecimal {
    assertEquals(mat.n, 2, "")
    val mc = MathContext(mat.expeps)
    val b = -mat[0, 0] - mat[1, 1]
    val c = mat[0, 0].multiply(mat[1, 1], mc) - mat[0, 1].multiply(mat[1, 0], mc)
    val d = b.pow(2, mc) - c.multiply(BigDecimal(4), mc)
    if (d < BigDecimal(0)) return BigDecimal(0)
    if (d == BigDecimal(0)) return b.divide(BigDecimal(-2), mc)
    val x1 = (b + d.sqrt(mc)).divide(BigDecimal(-2), mc)
    val x2 = (b - d.sqrt(mc)).divide(BigDecimal(-2), mc)
    if ((x1 - mat[1, 1]).abs() > (x2 - mat[1, 1]).abs())
        return x2
    return x1
}

fun QR_algo_free_diag_faster2(matrix: thee_diag_matrix, eps: Double,Q_out: Matrix): Pair<thee_diag_matrix, Matrix> {
    var E = thee_diag_matrix(
        matrix.n,
        Vector(Array(matrix.n - 1) { BigDecimal(0) }),
        Vector(Array(matrix.n) { BigDecimal(1) }),
        Vector(Array(matrix.n - 1) { BigDecimal(0) })
    )
    if (matrix.n < 2) return QR_algo_free_diag(matrix, eps)
    val s = take_s(sub_matrix(zdvig_m = matrix.n - 2, zdvig_n = matrix.n - 2, v_matrix = matrix))
    matrix -= E * s
    var (Q, R) = QR_razl_G_free_diag2(matrix,Q_out)
    matrix += E * s
//    var Q_out = (Q as my_Matrix).copy()
//    Q_out.copy_From(Q.transpose())
    var A = thee_diag_matrix(R.n)
    A.copy_From(R.mult_r_for_3_diag(Q) + (E3(R.n) * s))
    while (sub_matrix(zdvig_n = matrix.n - 1, zdvig_m2 = 1, v_matrix = A).norm() > BigDecimal(eps)) {
        val s = take_s(sub_matrix(zdvig_m = matrix.m - 2, zdvig_n = matrix.n - 2, v_matrix = matrix))
        A.copy_From(A - E * s)
        val (Q, R) = QR_razl_G_free_diag2(A,Q_out)
        A.copy_From(R * Q + E * s)
//        Q_out = Q_out * Q
//        Q_out.copy_From(Q_out * Q)

    }
    val m_in = thee_diag_matrix(A.n - 1)
    m_in.copy_From(sub_matrix(zdvig_n2 = 1, zdvig_m2 = 1, v_matrix = A))
    val (D2, Q2) = QR_algo_free_diag_faster2(m_in, eps,Q_out)
    val A_ref = sub_matrix(zdvig_n2 = 1, zdvig_m2 = 1, v_matrix = A)
    A_ref.copy_From(D2)

    //    Q_out = Q_out * Q2
//    val Q_out_ref = sub_matrix(zdvig_n2 = 1, zdvig_m2 = 1, v_matrix = Q_out)
//    Q_out_ref.copy_From(Q_out_ref * Q2)
//    Q_out.copy_From(Q_out * Q2)
    return Pair(A, Q_out)
}

fun QR_algo_free_diag_faster(matrix: thee_diag_matrix, eps: Double): Pair<thee_diag_matrix, Matrix> {
    val (D,Q)=  QR_algo_free_diag_faster2(matrix,eps,E(matrix.n,matrix.m))
    return Pair(D,Q.transpose())
}

fun super_faste_diag(matrix: Matrix,eps: Double):thee_diag_matrix{
    val (A,_) = free_diag(matrix)
    val (A2,_) = QR_algo_free_diag_faster(A,eps)
    return A2
}