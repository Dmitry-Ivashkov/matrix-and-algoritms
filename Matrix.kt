import java.lang.Integer.max
import java.math.BigDecimal
import java.math.MathContext
import kotlin.math.cos
import kotlin.math.min
import kotlin.math.sin
import kotlin.test.assertEquals

interface Matrix {  // транспонированное хранение данных
    val n: Int
        get
    val m: Int
        get
    val expeps
        get() = 10

    operator fun get(i: Int, j: Int): BigDecimal
    operator fun set(i: Int, j: Int, value: BigDecimal)
    fun gett(i: Int): IVector
    fun sett(i: Int, value: Vector)

//    override fun toString(): String {
//        return buildString(6 * m * n) {
//            for (i in 0 until n) {
//                for (j in 0 until m) {
//                    this.append(get(i, j).setScale(3, BigDecimal.ROUND_HALF_UP).toString())
//                    this.append(' ')
//                }
//                this.append('\n')
//            }
//            this.append('\n')
//        }
//    }

    operator fun times(matrix: Matrix): Matrix {
        assertEquals(m, matrix.n, "n!=m v *")
        val new = my_Matrix(this.n, matrix.m)
        for (i in 0 until this.n) {
            for (j in 0 until matrix.m) {
                for (k in 0 until this.m) {
                    new.set(i, j, new.get(i, j) + this.get(i, k).multiply(matrix.get(k, j), MathContext(expeps)))
                }
            }
        }
        return new
    }

    operator fun times(matrix: thee_diag_matrix): Matrix {
        assertEquals(m, matrix.n, "n!=m v *")
        val new = my_Matrix(this.n, matrix.m)
        for (i in 0 until this.n) {
            val k_beg: Int = max(0, i - 1)
            val k_end: Int = min(i + 2, this.m)
            for (j in 0 until this.m) {
                for (k in k_beg.until(k_end)) {
                    new[i, j]+= this[i, k].multiply(matrix[k, j], MathContext(expeps))
                }
            }
        }
        return new
    }

    operator fun timesAssign(s: BigDecimal) {
        for (i in 0 until this.n) {
            for (j in 0 until this.m) {
                set(i, j, get(i, j).multiply(s, MathContext(expeps)))
            }
        }
    }

    operator fun plusAssign(matrix: Matrix) {
        assertEquals(Pair(n, m), Pair(matrix.n, matrix.m), "n!=n v +=")

        for (i in 0 until this.n) {
            for (j in 0 until this.m) {
                this[i, j] = get(i, j) + matrix.get(i, j)
            }
        }
    }

    operator fun plus(matrix: Matrix): Matrix {
        assertEquals(Pair(n, m), Pair(matrix.n, matrix.m), "n!=n v +")
        val n = my_Matrix(this.n, this.m)
        for (i in 0 until this.n) {
            for (j in 0 until this.m) {
                n.set(i, j, get(i, j) + matrix.get(i, j))
            }
        }
        return n
    }

    operator fun minusAssign(matrix: Matrix) {
        assertEquals(Pair(n, m), Pair(matrix.n, matrix.m), "n!=n v +=")

        for (i in 0 until this.n) {
            for (j in 0 until this.m) {
                this[i, j] = get(i, j) - matrix.get(i, j)
            }
        }
    }

    operator fun minus(matrix: Matrix): Matrix {
        assertEquals(Pair(n, m), Pair(matrix.n, matrix.m), "n!=n v +")
        val n = my_Matrix(this.n, this.m)
        for (i in 0 until this.n) {
            for (j in 0 until this.m) {
                n.set(i, j, get(i, j) - matrix.get(i, j))
            }
        }
        return n
    }

    operator fun times(vector: Vector): Vector {
        assertEquals(m, vector.n, "m!=n v *vector")

        val new = Vector(this.n)
        for (i in 0 until this.n) {
            for (j in 0 until vector.m) {
                for (k in 0 until this.m) {
                    new.set(i, j, new.get(i, j) + this.get(i, k).multiply(vector.get(k, j), MathContext(expeps)))
                }
            }
        }
        return new
    }

    fun transpose(): Matrix {
        var new = my_Matrix(this.m, this.n)
        for (i in 0 until this.n) {
            for (j in 0 until this.m) {
                new.set(j, i, this.get(i, j))
            }
        }
        return new
    }

    fun norm(): BigDecimal {
        var nor = BigDecimal(0)
        for (i in 0 until this.n) {
            for (j in 0 until this.m) {
                nor += get(i, j).pow(2, MathContext(expeps))
            }
        }
        nor = nor.sqrt(MathContext(expeps))
        return nor
    }

    fun norm_ger(): BigDecimal {
        var nor = BigDecimal(0)
        for (i in 0 until this.n) {
            var v = BigDecimal(0)
            for (j in 0 until this.m) {
                if (i != j)
                    v += get(i, j).abs()
            }
            nor = nor.max(v)
        }
        return nor
    }

    open operator fun times(s: BigDecimal): Matrix {
        val n = my_Matrix(this.n, this.m)
        for (i in 0 until this.n) {
            for (j in 0 until this.m) {
                n.set(i, j, get(i, j).multiply(s, MathContext(expeps)))
            }
        }
        return n
    }

    operator fun times(h: H): Matrix {  //fast
        return this - ((this * (h.vector * BigDecimal(2))) * (h.vector.transpose()))
    }

    operator fun times(h: H2): Matrix {  //fast
        return ((this * (h.vector * BigDecimal(-2))) * (h.vector.transpose()))
    }

    fun mult_l(g: G) {  //fast
        for (i in 0 until m) {
            var v_i_ = get(g.i_, i)
            var v_j_ = get(g.j_, i)
            set(g.i_, i, g.c.multiply(v_i_, MathContext(expeps)) + g.s.multiply(v_j_, MathContext(expeps)))
            set(g.j_, i, g.c.multiply(v_j_, MathContext(expeps)) - g.s.multiply(v_i_, MathContext(expeps)))
        }
    }
    fun mult_r_for_3_diag(matrix: Matrix) : thee_diag_matrix{  //fast
        assertEquals(m, matrix.n, "n!=m v mult_r_for_3_diag")
        val new = thee_diag_matrix(this.n)
        for (i in 0 until this.n) {
            val k_beg: Int = max(0, i - 1)
            val k_end: Int = min(i + 2, this.m)
            for (k in k_beg.until(k_end)) {
                for (j in 0 until matrix.m) {
                    new[i, k] += this[i, j].multiply(matrix[j, k], MathContext(expeps))
                }
            }
        }
        return new
    }

    fun mult_l(h: H) {  //fast
        this -= (h.vector * BigDecimal(2)) * (h.vector.transpose() * this)
    }

    fun mult_lr(h: H2) {  //fast
        val matrix = h * this
        val matrix1 = this * h
        val matrix2 = matrix * h
        this += matrix + matrix1 + matrix2
    }

    fun copy(): Matrix {
        return this * BigDecimal(1)
    }

    fun copy_From(matrix: Matrix) {
        assertEquals(n, matrix.n, "n!=n v copy_From")
        for (i in 0 until this.n) {
            for (j in 0 until this.m) {
                this[i, j] = matrix[i, j]
            }
        }
    }

    operator fun times(cv: complex_vector): complex_vector {
        return complex_vector(this * cv.re, this * cv.im)
    }

    fun norm_with_matrix(cv: complex_vector): complex_namber {
        return complex_namber(
            (cv.re.transpose() * (this * cv.re))[0] + (cv.im.transpose() * (this * cv.im))[0],
            (cv.re.transpose() * (this * cv.im))[0] - (cv.im.transpose() * (this * cv.re))[0]
        )
    }
}

open class my_Matrix(var data: Array<Array<BigDecimal>>) : Matrix {  // транспонированное хранение данных
    override val n: Int = data.elementAtOrNull(0)?.size ?: 0
    override val m: Int = data.size
    override val expeps = 10

    constructor(n: Int = 0, m: Int = 0) : this(Array<Array<BigDecimal>>(m) { Array<BigDecimal>(n) { BigDecimal(0.0) } })
    constructor(arrays: Array<Array<Double>>) : this(arrays.map { array ->
        array.map { elem -> BigDecimal(elem) }.toTypedArray()
    }.toTypedArray())

    override fun toString(): String {
        return buildString(6 * m * n) {
            for (i in 0 until n) {
                for (j in 0 until m) {
                    this.append(get(i, j).setScale(3, BigDecimal.ROUND_HALF_UP).toString())
                    this.append(' ')
                }
                this.append('\n')
            }
            this.append('\n')
        }
    }

    override operator fun get(i: Int, j: Int): BigDecimal {
        return data[j][i]
    }

    override fun gett(i: Int): Vector {
        return Vector(data[i])
    }

    override operator fun set(i: Int, j: Int, bigDecimal: BigDecimal) {
        data[j][i] = bigDecimal
    }

    override fun sett(i: Int, vector: Vector) {
        data[i] = vector.data[0]
    }

}

open class E(n: Int, m: Int) : my_Matrix(n, m) {
    init {
        for (i in 0 until (min(n, m))) {
            set(i, i, BigDecimal(1))
        }
    }
}
open class E3(n: Int) : thee_diag_matrix(n, Vector(Array(n-1) { 0.0 }), Vector(Array(n) { 1.0 }), Vector(Array(n-1) { 0.0 })) {}

interface IVector : Matrix {

    operator fun get(n: Int): BigDecimal

    operator fun set(n: Int, bigDecimal: BigDecimal)

    override operator fun times(s: BigDecimal): Vector {
        val n = Vector(this.n)
        for (i in 0 until this.n) {
//            n.set(i, get(i).multiply(s, MathContext(expeps)))
            n[i] = get(i).multiply(s, MathContext(expeps));
        }
        return n
    }

    operator fun plus(vector: Vector): Vector {
        assertEquals(n, vector.n, "n!=n v +v")
        val n = Vector(this.n)
        for (i in 0 until this.n) {
            n.set(i, get(i) + vector.get(i))
        }
        return n
    }

    operator fun minus(vector: Vector): Vector {
        assertEquals(n, vector.n, "n!=n v -v")

        val n = Vector(this.n)
        for (i in 0 until this.n) {
            n.set(i, get(i) - vector.get(i))
        }
        return n
    }

    override fun copy(): IVector {
        return this * BigDecimal(1)
    }
}

class Vector(data: Array<BigDecimal>) : my_Matrix(Array(1) { data }), IVector {
    constructor(n: Int = 0) : this(Array<BigDecimal>(n) { BigDecimal(0.0) })
    constructor(array: Array<Double>) : this(array.map { elem -> BigDecimal(elem) }.toTypedArray()) {}

    override operator fun get(n: Int): BigDecimal {
        return data[0][n]
    }

    override operator fun set(n: Int, bigDecimal: BigDecimal) {
        data[0][n] = bigDecimal
    }

}

class L(n: Int, m: Int) : my_Matrix(n, m) { //  нижне треуг
    fun obr(vector: Vector): Vector {
        val new_v = Vector(m)
        for (i in 0 until m) {
            var z = vector[i]
            for (j in 0 until i) {
                z -= get(i, j).multiply(new_v[j], MathContext(expeps))
            }
//            z/=get(i,i)
            z = z.divide(get(i, i), MathContext(expeps))
            new_v[i] = z
        }
        return new_v
    }
}
class U(n: Int, m: Int) : my_Matrix(n, m) { //  верхне треуг
    fun obr(vector: Vector): Vector {
        val new_v = Vector(m)
        for (i in (m-1) downTo  0) {
//            print(i)
            var z = vector[i]
            for (j in i until n) {
                z -= get(i, j).multiply(new_v[j], MathContext(expeps))
            }
//            z/=get(i,i)
            z = z.divide(get(i, i), MathContext(expeps))
            new_v[i] = z
        }
        return new_v
    }
}

class G(n: Int = 0, m: Int = 0, val i_: Int, val j_: Int, val c: BigDecimal, val s: BigDecimal) : E(n, m) {

    constructor(n: Int = 0, m: Int = 0, i_: Int, j_: Int, phi: Double) : this(
        n,
        m,
        i_,
        j_,
        BigDecimal(cos(phi)),
        BigDecimal(sin(phi))
    )

    init {
        if ((n != 0) && (m != 0)) {
            set(i_, i_, c)
            set(j_, j_, c)
            set(i_, j_, s)
            set(j_, i_, -s)
        }
    }
}

class H(n: Int = 0, val vector: Vector) : E(n, n) {

    constructor(vector: Vector) : this(0, vector)

    init {
        if ((n != 0) && (m != 0)) {
            this -= ((vector * BigDecimal(2)) * vector.transpose())
        }
    }

    override operator fun times(matrix: Matrix): Matrix {  //fast
        return matrix - (vector * BigDecimal(2)) * (vector.transpose() * matrix)
    }

}

class H2(n: Int = 0, val vector: Vector) : E(n, n) {

    constructor(vector: Vector) : this(0, vector)

    init {
        if ((n != 0) && (m != 0)) {
            this -= ((vector * BigDecimal(2)) * vector.transpose())
        }
    }

    override operator fun times(matrix: Matrix): Matrix {  //fast
        return (vector * BigDecimal(-2)) * (vector.transpose() * matrix)
    }

}

class sub_matrix(
    val zdvig_n: Int = 0,
    val zdvig_m: Int = 0,
    val zdvig_n2: Int = 0,
    val zdvig_m2: Int = 0,
    val v_matrix: Matrix
) :
    Matrix {
    override val m: Int
        get() = v_matrix.m - zdvig_m - zdvig_m2
    override val n: Int
        get() = v_matrix.n - zdvig_n - zdvig_n2

    override fun get(i: Int, j: Int): BigDecimal {
        return v_matrix[i + zdvig_n, j + zdvig_m]
    }

    override fun set(i: Int, j: Int, value: BigDecimal) {
        v_matrix[i + zdvig_n, j + zdvig_m] = value
    }

    override fun gett(i: Int): IVector {
        return sub_vector(zdvig_n, 0, zdvig_n2, 0, v_matrix.gett(i + zdvig_m))

    }

    override fun sett(i: Int, vector: Vector) {
        v_matrix.sett(i + zdvig_m, vector)
    }

    override fun toString(): String {
        return buildString(6 * m * n) {
            for (i in 0 until n) {
                for (j in 0 until m) {
                    this.append(get(i, j).setScale(3, BigDecimal.ROUND_HALF_UP).toString())
                    this.append(' ')
                }
                this.append('\n')
            }
            this.append('\n')
        }
    }
}

class sub_vector(
    val zdvig_n: Int = 0,
    val zdvig_m: Int = 0,
    val zdvig_n2: Int = 0,
    val zdvig_m2: Int = 0,
    val v_vector: IVector
) :
    IVector {
    override val m: Int
        get() = v_vector.m - zdvig_m - zdvig_m2
    override val n: Int
        get() = v_vector.n - zdvig_n - zdvig_n2

    override fun get(i: Int, j: Int): BigDecimal {
        return v_vector[i + zdvig_n, j + zdvig_m]
    }

    override fun set(i: Int, j: Int, value: BigDecimal) {
        v_vector[i + zdvig_n, j + zdvig_m] = value
    }

    override fun get(i: Int): BigDecimal {
        return v_vector[i + zdvig_n]
    }

    override fun set(i: Int, value: BigDecimal) {
        v_vector[i + zdvig_n] = value
    }

    override fun gett(i: Int): IVector {
        return v_vector.gett(i + zdvig_m)
    }

    override fun sett(i: Int, vector: Vector) {
        v_vector.sett(i + zdvig_m, vector)
    }

    override fun toString(): String {
        return buildString(6 * m * n) {
            for (i in 0 until n) {
                for (j in 0 until m) {
                    this.append(get(i, j).setScale(3, BigDecimal.ROUND_HALF_UP).toString())
                    this.append(' ')
                }
                this.append('\n')
            }
            this.append('\n')
        }
    }
}

open class thee_diag_matrix(override val n: Int, var vector_r: Vector, var vector_m: Vector, var vector_l: Vector) :
    Matrix {
    override val m: Int
        get() = n

    constructor(n: Int) : this(
        n, Vector(Array(n - 1) { BigDecimal(0) }),
        Vector(Array(n) { BigDecimal(0) }),
        Vector(Array(n - 1) { BigDecimal(0) })
    )

    override fun get(i: Int, j: Int): BigDecimal {
        if (i == j) return vector_m[i]
        if (i == j + 1) return vector_l[j]
        if (i == j - 1) return vector_r[i]
        return BigDecimal(0)
    }

    override fun set(i: Int, j: Int, value: BigDecimal) {
        if (i == j) vector_m[i] = value
        if (i == j + 1) vector_l[j] = value
        if (i == j - 1) vector_r[i] = value
    }

    override fun gett(i: Int): Vector {
        throw Exception("don't use")
        return vector_m
    }

    override fun sett(i: Int, vector: Vector) {
        throw Exception("don't use")
    }

    override fun toString(): String {
        return buildString(6 * m * n) {
            for (i in 0 until n) {
                for (j in 0 until m) {
                    this.append(get(i, j).setScale(3, BigDecimal.ROUND_HALF_UP).toString())
                    this.append(' ')
                }
                this.append('\n')
            }
            this.append('\n')
        }
    }

    override operator fun times(matrix: Matrix): Matrix {
        assertEquals(m, matrix.n, "n!=m v *")
        val new = my_Matrix(this.n, matrix.m)
        for (i in 0 until this.n) {
            val k_beg: Int = max(0, i - 1)
            val k_end: Int = min(i + 2, this.m)
            for (j in k_beg.until(k_end)) {
                for (k in 0 until matrix.m) {
                    new[i, k] += this[i, j].multiply(matrix[j, k], MathContext(expeps))
                }
            }
        }
        return new
    }

//    override operator fun times(matrix: thee_diag_matrix): Matrix {
//        assertEquals(m, matrix.n, "n!=m v *")
//        val new = thee_diag_matrix(n)
//        for (i in 0 until this.n) {
//            val k_beg: Int = max(0, i - 1)
//            val k_end: Int = min(i + 2, this.m)
//            for (j in k_beg.until(k_end)) {
//                for (k in 0 until matrix.m) {
//                    new[i, k] += this[i, j].multiply(matrix[j, k], MathContext(expeps))
//                }
//            }
//        }
//        return new
//    }

    override fun copy_From(matrix: Matrix) {
        for (i in 0 until this.n) {
            val k_beg: Int = max(0, i - 1)
            val k_end: Int = min(i + 2, this.m)
            for (k in k_beg.until(k_end)) {
                this[i, k] = matrix[i, k]
            }
        }
    }

    override fun copy(): thee_diag_matrix {
        return thee_diag_matrix(n, vector_r.copy() as Vector, vector_m.copy() as Vector, vector_l.copy() as Vector)
    }

    override fun transpose(): thee_diag_matrix {
        val new = this.copy()
        val tmp = new.vector_l
        new.vector_l = new.vector_r
        new.vector_r = tmp
        return new
    }
}


open class complex_vector(val re: Vector, val im: Vector) {

    fun norm(): BigDecimal {
        var nor = re.norm().pow(2, MathContext(re.expeps)) + im.norm().pow(2, MathContext(re.expeps))
        return nor.sqrt(MathContext(re.expeps))
    }

    operator fun times(bigDecimal: BigDecimal): complex_vector {
        return complex_vector(re * bigDecimal, im * bigDecimal)
    }

    operator fun times(cn: complex_namber): complex_vector {
        return complex_vector(re * cn.re - im * cn.im_, im * cn.re_ + re * cn.im_)
    }

    operator fun plus(cv: complex_vector): complex_vector {
        return complex_vector(re + cv.re, im + cv.im)
    }

    operator fun minus(cv: complex_vector): complex_vector {
        return complex_vector(re - cv.re, im - cv.im)
    }

    operator fun timesAssign(bigDecimal: BigDecimal) {
        re *= bigDecimal
        im *= bigDecimal
    }

    override fun toString(): String {
        return buildString(6 * re.n) {
            for (i in 0 until re.n) {
                this.append("(re:")
                this.append(re[i].setScale(3, BigDecimal.ROUND_HALF_UP).toString())
                this.append(" im:")
//                this.append('(')
//                this.append(re[i].setScale(3, BigDecimal.ROUND_HALF_UP).toString())
//                this.append(' ')
                this.append(im[i].setScale(3, BigDecimal.ROUND_HALF_UP).toString())
                this.append(')')

                this.append('\n')
            }
            this.append('\n')
        }
    }
}

class complex_namber(re_: BigDecimal, im_: BigDecimal) : complex_vector(Vector(arrayOf(re_)), Vector(arrayOf(im_))) {
    var im_: BigDecimal
        get() = im[0]
        set(value) {
            im[0] = value
        }
    var re_: BigDecimal
        get() = re[0]
        set(value) {
            re[0] = value
        }
}

fun myMod(a: Int, n: Int): Int {
    return (a % n + n) % n
}

class z13_1(n_: Int) : my_Matrix(n_ * n_, n_ * n_) {
    init {
        for (i in 0 until this.n) {
            for (j in 0 until this.m) {
                val x1=i%n_
                val x2=j%n_
                val y1=i/n_
                val y2=j/n_
                if ((x1==x2)&&((y1==myMod(y2+2*x2,n_)||y1==myMod(y2-2*x2,n_)||y1==myMod(y2+2*x2+1,n_)||y1==myMod(y2-2*x2-1,n_))))
                    this[i, j] = BigDecimal(1)
                if ((y1==y2)&&((x1==myMod(x2+2*y2,n_)||x1==myMod(x2-2*y2,n_)||x1==myMod(x2+2*y2+1,n_)||x1==myMod(x2-2*y2-1,n_))))
                    this[i, j] = BigDecimal(1)

            }
        }
    }
}
class z13_2(val n_: Int) : my_Matrix(n_+1, n_ +1) {
    init {
        for (i in 0 until this.n-1) {
            for (j in 0 until this.m-1) {
                if ((i==j-1)||(i==j+1)||(myMod(i*j, n_)==1))
                    this[i, j] = BigDecimal(1)
            }
        }
        this[0,n_]=BigDecimal(1)
        this[n_,0]=BigDecimal(1)
        this[n_,n_]=BigDecimal(2)
    }
}
