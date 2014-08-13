# Intro to Linear Algebra
## Justin Lynd / Richard Voepel
## Spring / Summer 2014

---
## Table of Contents

1. [Matrices and Vectors](#anchor1.1)
	1. [Linear Combinations](#anchor1.2)
	- [Gaussian Elimination](#anchor1.4)
		- [Rank and Nullity](#anchor1.41)
	- [Linear Independence](#anchor1.7)
2. [Invertibility](#anchor2.3)
	- [Inverses](#anchor2.4)
3. [Determinants]($anchor3.1)
4. [Subspaces](#anchor4.1)
	- [Bases](#anchor4.2)
5. [Eigen-things](#anchor5)
	1. [Characterstic Polynomial](#anchor5.2)
	- [Diagonalization of Matrices](#anchor5.3)
6. [Orthogonality](#anchor6)



---
### 1/22/14

#### Concrete
- Study of solutions of/how to solve systems of linear equations (Sect. 1.1, 2.6, 3)

#### Abstract
- (Sect. 4 & 5)

#### Usefulness
- Can be used to "linearize" mathematical objects to make it easier to predict behavior
- Ex:
	\\[\begin{aligned}
	 x-y &=3 \\\
	 2x+3y &=1 \\\
	\end{aligned} \\]
	
- Can be solved with either substitution or elimination
- **3-dimensions**:
	\\[\begin{aligned}
	x+2y-z &=1 \\\
	2x+y+z &=0 \\\
	\end{aligned} \\]
- **5-dimensions**
	$$2x_1 + 3x_2 - x_3 + x_4 - x_5$$
	
#### From Systems of Equations to Matrices and Vectors

Ex:
	\\[\begin{aligned}
	&2x_1 + x_2 + x_3 &= 1 \\\
	&4x_1 + x_2  &= -2 \\\
	-&2x_1 + 2x_2 + x_3 &= 7 \\\
	\end{aligned} \\]
	
Unknowns:
	$$\vec{x} = 
	\begin{bmatrix} x_1 \\\ x_2 \\\ x_3 \\\ \end{bmatrix}$$

	
#### Dot Product of two vectors

$$\vec{u} \bullet \vec{v} = \begin{bmatrix} u_1 \end{bmatrix}$$

##### Definition:
Given an *m x n* matrix $$$A =[\vec a_1 \dots \vec a_n]$$$
 and a *n x 1* vector $$$\vec{x} = \begin{bmatrix} x_1 \\\ x_2 \\\ \vdots \\\ x_n \end{bmatrix}$$$ the {matrix $$$\bullet$$$ vector} product is the sum scalar
 $$ x_1\vec a_1 = x_2\vec a_2 +\dotso+x_n\vec a_n$$
 
 ---
##  Text. Sect. 1.1

### [Matrices and Vectors](id:anchor1.1)

#### Definition:
A **matrix** is a rectangular array of scalars. If the matrix has *m* rows and *n* columns, we say that the *size* of the matrix is **m by n**, written *m x n*. The matrix is **square** if $$$m = n$$$. The scalar in the *ith* row and *jth* column is called the $$$(i,j)\mathrm{-entry}$$$ of the matrix.

- The **scalar multiple** $$$cA$$$ is the *m x n* matrix whose entries are *c* times the corresponding entries of *A*; $$$cA$$$ is the *m x n* matrix whose $$$(i,j)$$$-entry is $$$ca_{ij}$$$. (-1)A is -A and 0A is denoted by *O*, the **zero matrix**, where each entry is 0.

#### Properties of Matrix Addition and Scalar Multiplication (Theorem 1.1):
Let A, B, and C be *m x n* matrices, and let *s* and *t* be any scalars. Then:

1. $$$A+B = B+A$$$
2. $$$A+(B+C)=(A+B)+C$$$
3. $$$A+0\text{ [zero matrix]} = A$$$
4. $$$A + (-A) = 0$$$
5. $$$(st)A = s(tA)$$$
7. $$$s(A+B) = sA + sB$$$
8. $$$(s+t)A = sA + tA$$$

Proof of 7:
1. $$$[(s+t)A]\_{ij} = sA\_{ij} + tA\_{ij} = [sA]\_{ij} + [tA]\_{ij}$$$

 ---
 
###  1/27/14

### Transpose
- The **transpose** of an *m x n* matrix *A* is the *n x m* matrix denoted by $$$A^T$$$ whose $$$(i,j)$$$-entry is the $$$(j,i)$$$-entry of *A*:
\\[\begin{align}
&A=(a\_{ij}) \\\
&A^T=(a_{ji}) \\\
\\end{align} \\]
- If we take a real number to be a $$$1\times 1$$$ matrix, then its transpose is itself

#### Properties of Transpose (Thm 1.2):

Suppose A, B are $$$n\times m$$$ matrices, and $$$s\in\mathbb{R}$$$

1. $$$(A+B)^T = A^T + B^T$$$
2. $$$(sA)^T = sA^T$$$
3. $$$(A^T)^T = A$$$

	
### Vectors

- A matrix $$$\mathbf u$$$ that has exactly one row is a **row vector**; likewise with columns.
	- Usually be dealing with **column vectors**.
- Vectors are denoted with boldface; **u** and **v** or $$$\vec u$$$ and $$$\vec v$$$
- The set of all column vectors with *n* components is written as $$$\mathbb R^n$$$
- Vector entries are called **components**
- Vectors have both *magnitude and direction*
- 
- You can add and multiply vectors by scalars.
	- These operations are **vector addition** and *scalar multiplication*
- **0** is the **zero vector**, which obeys the following rules:
\\[ \begin{align}
\mathbf u + \mathbf 0 &= \mathbf u \\\
0\mathbf u &= \mathbf 0 \\\
\forall \mathbf u &\in \mathbb R^ n
\end{align} \\]

- For any *m x n* matrix $$$A$$$, its *jth* column is $$\mathbf a_ j = \begin{bmatrix} a\_{1j} \\\  a\_{2j} \\\ \vdots \\\  a\_{mj} \end{bmatrix}$$
- Vectors can also be represented geometrically as directed line segments/arrows
	- If $$$\mathbf v = \begin{bmatrix} a \\\ b \\\ \end{bmatrix}$$$ is a vector in $$$\mathbb R^2$$$ you can represent it as an arrow from the origin to the point $$$(a,b)$$$ in the $$$xy$$$-plane [insert graph thing in your mind]

### Linear Combinations

- **Definition:** Let $$$\vec u_1, \dotsc , \vec u_n \;\;\text{in}\;\; \mathbb{R}^m$$$ (set of all $$$m\times 1$$$ matrices with real number entries)
	- A *linear combination* of $$${\vec u_1, \dotsc ,\vec u_n}$$$ is a vector $$$\vec u \;\;\text{in}\;\;\mathbb{R}^m$$$ which is of the form
	$$\vec u = c_1\vec u_1 + c_2\vec u_2 + \dotso + c_n\vec u_n$$
	where $$$c_i$$$ is a scalar.
	
#### Examples
1. $$$\vec O$$$ is a linear combination of any collection of vectors $$\vec O = o\*\vec u_1 + \dotso + o*\vec u_n$$
2. $$$\vec u_1 = 1\*\vec u_2 + o*\vec u_2 + \dotso + o\*\vec u_n$$$ <br><br>
3. Is $$$\begin{bmatrix} 1 \\\ 2 \\\ 3\\\ \end{bmatrix}$$$ a linear comb. of $$$\begin{bmatrix} 1 \\\ 0 \\\ 3 \\\ \end{bmatrix}\;\;\text{and}\;\;\begin{bmatrix} -1 \\\ 1 \\\ 3 \\\ \end{bmatrix}$$$?<br><br>
	- i.e. are there weights $$$x_1, x_2,$$$ where $$x_1\begin{bmatrix} 1\\\0\\\3\\\ \end{bmatrix} + x_2 \begin{bmatrix} -1 \\\ 1 \\\ 3 \\\ \end{bmatrix} = \begin{bmatrix} 1 \\\ 2 \\\ 3\\\ \end{bmatrix}$$
	- i.e. are there solutions to this:
	\\[\begin{aligned}
	x_1 - x_2 = 1 \\\
	0 + x_2 = 2 \\\
	3x_1 + 3x_2 = 3 \\\
	x_2 = 2 \\\
	x_1 = 3 \\\
	\end{aligned}\\]
	
	\\[\begin{aligned}
	\let \vec w_1 &= s\vec u + t\vec v \\\
	\vec w_2 &= p\vec u + q\vec v \\\
	\therefore\vec w_1 + \vec w_2 &= (s+p)\vec u + (t + q)\vec v \\\
	\end{aligned}\\]
 	
###  Standard Vectors

- The **standard vectors** of $$$\mathbb R^n$$$ are $$\begin{bmatrix}1 \\\ 0 \\\ \vdots \\\ 0 \end{bmatrix}, \begin{bmatrix}0 \\\ 1 \\\ \vdots \\\ 0\end{bmatrix}, \begin{bmatrix} 0 \\\ 0 \\\ \vdots \\\ 1\end{bmatrix}$$

---
- One question asked in two different ways:
	- Does the system have a solution?
	- Is the RHS a linear combination of the *columns* of the coefficient matrix?
- **Definition:** $$$\mathbf{A} =[\vec a_1, \dotsc , \vec a_n]$$$,
	- Matrix-vector product 
	\\[\begin{aligned}
	&\begin{bmatrix} 1 & 2 & 3 \\\ 4 & 5 & 6 \\\ \end {bmatrix} \begin{bmatrix}7\\\8\\\9\\\ \end{bmatrix} = 7\begin{bmatrix}1\\\4\\\ \end{bmatrix} + 8\begin{bmatrix}2\\\5\\\ \end{bmatrix} + 9\begin{bmatrix}3\\\6\\\ \end{bmatrix} \\\
	&\mathbf{A} \bullet \vec x= x_1\vec a_1 + x_2\vec a_2 + \dotso + x_n\vec a_n \\\
	\end{aligned} \\]
	
#### Rotations
- What's $$$\mathbf {A_\theta}$$$?
	- An $$$n\times n$$$ matrix such that $$$\mathbf A_\theta\vec u = \vec v$$$
$$\mathbf A\_\theta = \begin{bmatrix} \cos\theta & -\sin\theta \\\ \sin\theta & \cos\theta\end{bmatrix}$$
- insert complicated graph and trig here

#### Thm 1.3
Suppose $$$A, B = m\times n$$$ matrices and $$$\vec u,\vec v\in \mathbb R^n$$$:

1. $$$A(\vec u + \vec v) = A\vec u + A\vec v$$$
2. $$$A(c\vec u) = cA\vec u$$$
3. $$$(A + B)\vec u = A\vec u + B\vec u$$$
4. $$$A\vec e_i = \vec a_i, A=[\vec a_1, \vec a_2,\cdots \vec a_n]$$$
5. $$$\text{If } A\vec w=B\vec w \;\forall \vec w \in\mathbb R^n \therefore A = B$$$


### Sect. 1.3

#### Elementary Row Ops/ Row Echelon Form

Ex:
\\[\begin{aligned}
&2x + 4y - 2z = 2 \\\
&4x + 9y - 3z = 8 \\\
-&2x - 3y + 7z = 10 \\\
\end{aligned} \\]

Ex:
\\[\begin{aligned}
x + 2y - z = 1 \\\
 y + z = 4 \\\
z = 2 \\\
\end{aligned} \\]

- Make into augmented matrix in *row echelon form*
- can be solved by "back substitution"
\begin{bmatrix}
2 & 4 & -2 &| &2 \\\
4 & 9 & -3 &| &8 \\\
-2 & -3 & 7 &| &10 \\\
\end{bmatrix}

\\[ \begin{aligned}
&r_1 + r_3 \to r_3 \\\
&r_2 + (-2)r_1 \to r_2 \\\
\end{aligned} \\]

\begin{bmatrix}
2 & 4 & -2 &| &2 \\\
0 & 1 & 1 &| &4 \\\
0 & 1 & 5 &| &12 \\\
\end{bmatrix}

---
###  Text. Sect. 1.2

## [Linear combinations](id:anchor1.2)

- A **linear combination** of vectors $$$\mathbf u_1, \mathbf u_2, \dotsc,\mathbf u_k$$$ is a vector of the form $$c_1\mathbf u_1 + c_2\mathbf u_2 + \dotso + c_k\mathbf u_k$$ where $$$c_1, c_2, \dotsc, c_k$$$ are scalars. These are called the **coefficients** of the linear combination.

### Matrix-Vector Products

- A **matrix-vector product** is the product of $$$m \times n$$$ matrix $$$A$$$ and $$$n \times 1$$$ vector **v**. $$$A\mathbf v$$$ is the linear combination of the columns of $$$A$$$ whose coefficients are the corresponding components of **v**. $$A\mathbf v = v_1\mathbf a_1+v_2\mathbf a_2+\dotso + v_n\mathbf a_n$$

- A matrix with nonnegative entries that sum to one is called a **stochastic matrix**

#### Identity Matrices

- For each positive integer $$$n$$$ the $$$n\times n\textbf{ identity matrix }I_n\text{ is the } n\times n$$$ matrix whose respective columns are the standard vectors $$$\mathbf e_1, \mathbf e_2, \dotsc,\mathbf e_n\text{ in } \mathbb R^n$$$

---
### Text. $$$\S$$$ 1.3

## Systems of Linear Equations

- A **linear equation** is an equation that can be written in the form $$a_1x_1 + a_2x_2 + \dotso + a_nx_n = b$$ where $$$a_1,a_2,\dotsc,a_n, b\in \mathbb R$$$.
- A **system of linear equations** is a set of *m* linear equations in the same *n* variables, where *m* and *n* are positive integers. This can be written in the following:

\\[\begin{aligned}
a\_{11}x\_1+a\_{12}x\_2+&\dotso+a\_{1n}x\_n = b\_1 \\\
a_{21}x_1 + a\_{22}x\_2+&\dotso+a\_{2n}x\_n = b\_2 \\\
&\vdots \\\
a\_{m1}x\_{2} + a\_{m2}x\_2+&\dotso + a\_{mn}x\_n = b\_m \\\
\end{aligned}\\]

- 
$$\mathrm{span}(v_1,v_2,\dotsc,v_n) = \\{c_1v_1 + c_2v_2 + \dotso c_nv_n\;|\;c_i\in\mathbb R\text{ for }1\leq i\leq n\\}$$

---
### 7/9/14

The **augmented matrix** associated with the equation $$$A\vec x=\vec b$$$, is the matrix $$$[\vec a_1, \vec a_2 \cdots \vec a_n \vec b]$$$

**Ex:**
\\[\begin{aligned}
x_1 - 2x_2 - x_3 = 3 \\\
3x_1 - 6x_2 - 5x_3 = 3 \\\
2x_1 - x_2 + x_3 = 0 \\\
\end{aligned}\\]
put into *augmented matrix* form:
\\[\begin{bmatrix}
1 & -2 & -1 & 3 \\\
3 & -6 & -5 & 3 \\\
2 & -1 & 1 & 0 \\\
\end{bmatrix}\\]

$$r_2-3r_1\to r_2$$

\begin{bmatrix}
1 & -2 & -1 & 3 \\\
0 & 0 & -2 & -6 \\\
2 & -1 & 1 & 0 \\\
\end{bmatrix}

$$r_3 - 2r_1 \to r_3$$

\begin{bmatrix}
1 & -2 & -1 & 3 \\\
0 & 0 & -2 & -6 \\\
0 & 3 & 3 & -6 \\\
\end{bmatrix}

$${1\over 3}r_3\to r_3$$

\begin{bmatrix}
1 & -2 & -1 & 3 \\\
0 & 0 & -2 & -6 \\\
0 & 1 & 1 & -2 \\\
\end{bmatrix}

$$ r_3 \leftrightarrow r_2$$

\begin{bmatrix}
1 & -2 & -1 & 3 \\\
0 & 1 & 1 & -2 \\\
0 & 0 & -2 & -6 \\\
\end{bmatrix}

$${r_3\over -2}\to r_3$$

\begin{bmatrix}
1 & -2 & -1 & 3 \\\
0 & 1 & 1 & -2 \\\
0 & 0 & 1 & 3 \\\
\end{bmatrix}

#### Thm 1.4
Every matrix can be transformed into one and only one matrix in rref by means of elementary row operations.

**Procedure:** Given a system of linear equations,:

1. Write the augmented matrix $$$\begin{bmatrix} A & \mathbf b \end{bmatrix}$$$
2. By iterating elementary row operations, ensure that subsequent augmented matrices adhere to criteria for rref
3. You will end up with a matrix $$$\begin{bmatrix}R & \mathbf c\end{bmatrix}$$$. If this matrix has a row with entries that are all 0 *except* for the last entry, the system of equations is *inconsistent*. Otherwirse, the system has at least one solution.

- Variables that are not *free* are **basic variables**

---
## $$$\S$$$ 1.4 - [Gaussian Elimination](id:anchor1.4)

Suppose that $$$R$$$ is the rref of a matrix $$$A$$$. The first nonzero entry in a row is the **leading entry** of that row. the $$$(i,j)$$$ position in $$$R$$$, containing one of the leading entries is called a **pivot position**. A column that contains a pivot position, is called a **pivot column**.

### Gaussian Elimination

Forwards pass:

1. Find pivot column, and row swap until you have a nonzero entry in that pivot position.
2. Zero out all entries in pivot column below pivot position.
3. Ignore previous row and repeat steps 1-2 for all subsequent rows.

Backwards pass:

1. Find bottom row, and zero out all entries above pivot position
2. Ignore previous row, and continue upwards repeating step 1 for all subsequent rows.

#### [Rank and Nullity](id:anchor1.41)
Let $$$R$$$ be the rref of matrix $$$A$$$. The **rank** of $$$A$$$ is the number of nonzero rows in $$$R$$$. The **nullity** of $$$A$$$ is defined to be the number of columns minus the *rank*.

1. The rank of $$$A$$$ = # of pivot columns in $$$R$$$
2. Nullity of $$$A$$$ = # of non-pivot columns in $$$R$$$.
	- Adding zero rows will not change *rank* or *nullity*
3. If $$$A$$$ has **full rank** (*rank* is equal to total number of columns), and is square, $$$\Rightarrow R = I_n$$$ 

- If $$$A\mathbf x = \mathbf b$$$ is consistent:
	- \# of basic variables = rank
	- \# free variables = nullity
	
Consider the following:
\begin{aligned}
&x_1 + x_2 + &x_3 &= 1 \\\
&x_1 + &3x_3 &= -2+s \\\
&x_1 - x_2 + &rx_3 &= 3 \\\
\end{aligned}

Rewrite as $$$A\mathbf x = \mathbf b$$$
\begin{bmatrix}
1 & 1 & 1 & 1 \\\
1 & 0 & 3 & (-2+s) \\\
1 & -1 & r & 3 \\\
\end{bmatrix}

Perform *Gaussian elimination*


$$r_2\leftrightarrow r_3$$

\begin{bmatrix}
1 & 1 & 1 & 1 \\\
1 & -1 & r & 3 \\\
1 & 0 & 3 & (-2+s) \\\
\end{bmatrix}

\begin{aligned}
r_2 - r_1\to r_2 \\\
r_3 - r_1\to r_3 \\\
\end{aligned}

\begin{bmatrix}
1 & 1 & 1 & 1 \\\
0 & -2 & (r-1) & 2 \\\
0 & -1 & 2 & (-3+s) \\\
\end{bmatrix}

$$r_2\leftrightarrow r_3$$

\begin{bmatrix}
1 & 1 & 1 & 1 \\\
0 & -1 & 2 & (-3+s) \\\
0 & -2 & (r-1) & 2 \\\
\end{bmatrix}

$$r_3 - 2r_2\to r_3$$
\begin{bmatrix}
1 & 1 & 1 & 1 \\\
0 & -1 & 2 & (-3+s) \\\
0 & 0 & (r-5) & 8-2s \\\
\end{bmatrix}

Three cases:

If $$$r=5,s=4$$$:
\begin{bmatrix}
1 & 1 & 1 & 1 \\\
0 & -1 & 2 & (-3+s) \\\
0 & 0 & 0 & 0 \\\
\end{bmatrix}
Infinitely many solutions, 3 unknowns for 2 equations.

If $$$r=5,s\neq 4$$$:
\begin{bmatrix}
1 & 1 & 1 & 1 \\\
0 & -1 & 2 & (-3+s) \\\
0 & 0 & 0 & (8-2s) \\\
\end{bmatrix}
No solution.

If $$$r\neq 5$$$:
\begin{bmatrix}
1 & 1 & 1 & 1 \\\
0 & -1 & 2 & (-3+s) \\\
0 & 0 & 1 & {8-2s\over r-5} \\\
\end{bmatrix}
$$\vdots$$
\begin{bmatrix}
1 & 0 & 0 & ? \\\
0 & 1 & 0 & ? \\\
0 & 0 & 1 & ? \\\
\end{bmatrix}
One unique solution.

#### Thm 1.5
The following are equivalent:

1. The equation $$$A\mathbf x = \mathbf b$$$ is consistent (has at least one solution)
2. The vector $$$\mathbf b$$$ is a linear combination of the columns of $$$A$$$.
3. The *rref* of $$$A$$$ has no row of the form $$$\begin{bmatrix}0 & 0 & 0 & \dots & 0 & d\end{bmatrix},\; d\neq 0$$$.

---
### 7/11/14

## $$$\S$$$ 1.7 - [Linear Independence](id:anchor1.7)

**Definition:** A matrix $$$A = \\{\vec u_1, \vec u_2, \dotsc , \vec u_k\\}$$$ is **linearly dependent** iff $$$A\mathbf x = \mathbf 0$$$

When determining linear independency, solve the matrix equation $$$A\mathbf x = \mathbf 0$$$.

Let $$$A$$$ be an $$$m\times n $$$ matrix. The following are equivalent:

1. The columns of A are *linearly independent*
2. The equation $$$A\mathbf x = \mathbf b$$$ has at *most* one solution for every $$$\mathbf b$$$.
3. The nullity of $$$A$$$ is 0.
4. The rank of $$$A$$$ is equal to the number of columns $$$n$$$.
5. The columns of rref of $$$A$$$ are distinct standard vectors in $$$\mathbb R^n$$$.
6. The only solution to $$$A\mathbf x = \mathbf 0$$$ is $$$\mathbf x = \mathbf 0$$$.
7. There is a pivot position in each column of $$$A$$$. 

The matrix equation $$$A\mathbf x = \mathbf 0$$$ is said to be a homogeneous equation.

$$$\[A \;\;\mathbf 0] = \begin{bmatrix} 1 & -4 & 2 & -1 & 2 & 0 \\\ 2 & -8 & 3 & 2 & -1 & 0 \end{bmatrix}$$$

End result:
\\[\begin{bmatrix}
4x_2 - 7x_4 + 8x_5 \\\
x_2 \\\
4x_4 - 5x_3 \\\
x_4 \\\
x_5 \\\
\end{bmatrix}\\]

$$= x_2\begin{bmatrix} 4 \\\ 1 \\\ 0 \\\ 0 \\\ 0 \end{bmatrix} + x_4\begin{bmatrix} -7 \\\ 0 \\\ 4 \\\ 1 \\\ 0 \end{bmatrix} + x_5 \begin{bmatrix} 8 \\\ 0 \\\ -5 \\\ 0 \\\ 1 \\\ \end{bmatrix}$$


**Fact:** The set of vectors produced by Gaussian Elimination are *linearly independent*

### Thm 1.8

Let $$$\\{\vec u_1, \dotsc,\vec u_k\\}\subseteq \mathbb R^n$$$. This set of vectors is *linearly dependent* iff $$$\vec u_i = \vec 0$$$, or $$$\\{\vec u_1,\dotsc,\vec u_j\\}$$$ is linearly dependent. Moreover, it will be the case that $$$\vec u_j = c_1\vec u_1 + \dotso + c\_{j-1} u\_{j-1}$$$.

1. A set consisting of a single vector is linearly independent only when it is a non-zero vector.
2. A set of two vectors \\{\vec u_1,\vec u_2\\} is linearly dependent iff $$$\vec u_1 = \vec 0$$$ or $$$\vec u_2 = \vec 0$$$ or if $$$\vec u_1 = c\vec u_2$$$ for some $$$c\in \mathbf R^2$$$.
3. Set $$$S=\\{\vec u_1,\dotsc,\vec u_k\\}$$$ be linearly independent, $$$\vec v\in \mathbb R^n$$$. If $$$\vec v \notin \mathrm{span}(S)$$$, then $$$\\{\vec u_1, \dotsc, \vec u_k, \vec v\\}$$$ is also linearly independent.
4. Every subset of $$$\mathbb R^n$$$ that contains *strictly* more than $$$n$$$ vectors, is linearly dependent.
5. If $$$S\subseteq \mathbb R^n$$$ and no vector can be removed from $$$S$$$ without altering the $$$\mathrm{span}(S)$$$, then $$$S$$$ is linearly independent.

---
### 7/14/14

## $$$\S$$$ 2.3 - [Invertibility](id:anchor2.3)

**Definition:** Let $$$A$$$ be an $$$n\times n$$$ matrix and suppose that $$$B$$$ is an $$$n\times n$$$ matrix that has the property $$$AB=I_n=BA$$$. Then we call $$$B$$$ the **inverse** of $$$A$$$.

$$$\begin{bmatrix} 1 & 0 & 0 \\\ 0 & 2 & 0 \\\ 0 & 0 & 3 \end{bmatrix}
\begin{bmatrix} 1 & 0 & 0 \\\ 0 & ½ & 0 \\\ 0 & 0 & ⅓ \end{bmatrix}
= \begin{bmatrix} 1 & 0 & 0 \\\ 0 & 1 & 0 \\\ 0 & 0 & 1 \end{bmatrix}$$$

\\[\begin{aligned}
&xr = s \\\
\implies &x = r^{-1} s \\\
&A\mathbf x = \mathbf b \\\
\implies &BA\mathbf x = B\mathbf b \\\
\implies &I_n \mathbf x = B\mathbf b \\\
\implies &\mathbf x = B\mathbf b \\\
\end{aligned}\\]

**Fact:** If $$$A$$$ has an inverse, $$$A\mathbf x = \mathbf b$$$ is consistent for all $$$\mathbf b$$$, and has exactly one solution.

#### Thm 2.1

The matrix inverse is unique.
**Proof:** Let $$$A$$$ be an $$$n\times n$$$ matrix, and assume that $$$B, C$$$ are both inverses for $$$A$$$

$$B=BI_n = B(AC) = (BA)C = I_n C = C$$

#### Thm 2.2

Let $$$A, B$$$ be $$$n\times n$$$ matrices. The following are true:

1. If $$$A^{-1}$$$ exists, $$$A^{-1}$$$ is invertible, and its inverse $$$(A^{-1})^{-1} = A$$$
2. If $$$A^{-1}, B^{-1}$$$ both exist, $$$(AB)^{-1} = B^{-1} A^{-1}$$$
	3. For invertible $$$n\times n$$$ matrcies $$$A_1,\dotsc, A_k$$$ $$(A_1 A_2\dotso A_k)^{-1} = A_k^{-1} A\_{k-1}^{-1}\dotso A_1^{-1}$$
3. $$$(A^{-1})^T = (A^T)^{-1}$$$, should $$$A^{-1}$$$ exist.

### Elementary Matrices

1. $$$P_{i,j}$$$ is the elementary matrix obtained by swapping rows $$$i$$$ and $$$j$$$ in the identity matrix.
2. $$$S_{r,i}$$$ is the elementary matrix obtained by scaling the $$$i^{th}$$$ row of the identity matrix by $$$r$$$.
3. $$$A_{i,j}$$$ is the elementary matrix obtained by placing an extra 1 in the $$$(i,j)^{th}$$$ position of the identity matrix, where $$$i\neq j$$$. This is equivalent to row $$$i$$$ plus row $$$j$$$ being stored in row $$$i$$$.

**Claim:** There exist $$$e_1, e_2, \dotso, e_k$$$ of elementary matrices, such that $$$e_k, e\_{k-1}, \dotso , e_1 A=R$$$  where $$$R$$$ is the rref of $$$A$$$. There is a *single* invertible matrix $$$P$$$ such that $$$PA=R$$$.

#### Column Correspondence Property

Given an $$$m\times n$$$ matrix $$$A$$$, Compute $$$R$$$.
If $$$\mathbf r_j$$$ is a linear combination of some other set of columns, $$$\\{\mathbf r_l\\}$$$, with scalars $$$\\{c_l\\}$$$, then $$$\mathbf a_j$$$ is a linear combination of $$$\\{\mathbf a_l\\}$$$ with scalars $$$\\{c_l\\}$$$.

#### Thm 2.4

The following statements are true:

1. the pivot columns of any matrix $$$A$$$ are linearly independent.
2. Each non-pivot columns of $$$A$$$ is a linear combination of the previous pivot columns of $$$A$$$, where the entries (in rref) are the scalars associated with the linear combination.

---
### 7/16/14

## $$$\S$$$ [2.4 - Inverses](id:anchor2.4)

#### Thm 2.5
Let $$$A$$$ be a $$$n\times n$$$ matrix. Then $$$A$$$ is invertible iff the rref of $$$A$$$ is $$$I_n$$$.

##### Algorithm for building an Inverse:

Before: $$$A\mathbf x = \mathbf b$$$, solve for $$$\mathbf x$$$.
Now: $$$AB=I_n$$$, solve for $$$B$$$:
	
- Consider $$$\begin{bmatrix} A & \mid & I \end{bmatrix}$$$
- Apply Gaussian elimination, find $$$\begin{bmatrix} R & \mid & B \end{bmatrix}$$$, where $$$R=\mathrm{rref}(A)$$$, and $$$B$$$ is some matrix.

### Thm 2.6

Let $$$A$$$ be an $$$n\times n$$$ matrix. The following are equivalent.

1. $$$A$$$ is invertible
2. $$$\mathrm{rref}(A)$$$ = $$$I_n$$$
3. rank$$$(A)=n$$$ (full rank)
4. span of the columns of $$$A$$$ is all of $$$\mathbb R^n$$$
5. $$$A\mathbf x = \mathbf b$$$ is consistent $$$\forall\mathbf b\in \mathbb R^n$$$
6. nullity$$$(A)=0$$$
7. There is some $$$n\times n$$$ matrix $$$B$$$, $$$\ni AB=I_n$$$
8. There is an $$$n\times n$$$ matrix $$$C\ni C=I_n$$$
9. The columns of $$$A$$$ are linearly independent
10. The only solution of $$$A\mathbf x = \mathbf 0$$$ is $$$\mathbf x = \mathbf 0$$$
11. $$$A=\\{E_k E\_{k-1} \dotso E_2 E_1\\}$$$ where each $$$E_j$$$ is an elementary matrix.

---
### 7/18/14
## $$$\S$$$ 2.5 - Block Multiplication

Suppose we wanted to calculate $$$A^2$$$ for $$A = \begin{bmatrix}
1 & 0 & 0 & 0 \\\
0 & 1 & 0 & 0 \\\
6 & 8 & 5 & 0 \\\
-7 & 9 & 0 & 5 \\\
\end{bmatrix}$$
If we wanted to do this with matrix multiplication it would be very tedious and cumbersome. However, we can partition the matrix into blocks and make the task much easier:

\\[\begin{align}
A &= \begin{bmatrix}
1 & 0 & \mid & 0 & 0 \\\
0 & 1 & \mid & 0 & 0 \\\
- & - & + & - & - \\\
6 & 8 & \mid & 5 & 0 \\\
-7 & 9 & \mid & 0 & 5 \\\
\end{bmatrix} 
= \begin{bmatrix}
I & \mid & 0 \\\
- & + & - \\\
B & \mid & 5I \\\
\end{bmatrix} \\\
\\\
A^2 &= \begin{bmatrix}
I & \mid & 0 \\\
- & + & - \\\
B & \mid & 5I \\\
\end{bmatrix}
\begin{bmatrix}
I & \mid & 0 \\\
- & + & - \\\
B & \mid & 5I \\\
\end{bmatrix} \\\
\\\
&= \begin{bmatrix}
I^2 & \mid & 0 \\\
- & + & - \\\
6B & \mid & (5I)^2 \\\
\end{bmatrix}
\end{align}\\]

## $$$\S$$$ 2.6 - LU Decomposition



---


### 3/12/12

## $$$\S$$$ 3.1 - [Determinants](id:anchor3.1)

\\[\begin{align}
&\det(A) = a\_{11}c\_{11} + a\_{12}c\_{12} +\dotso + a\_{1n}c\_{1n} \\\
&\text{where } c\_{ij} = (-1)^{i+j}\,\det\,(A\_{ij}) \\\
\end{align}\\]

If $$$\det(A)\neq 0,\; A$$$ is invertible.

#### Thm 3.2

- det(triangular matrix) = product of diag. entries

##### Notation:

$$\det(A) = \mid A\mid$$


So let's consider $$$\\{\mathbf u, \mathbf v\\} \leq \mathbb R^2$$$ that are linearly independent.
- Recall that **u** and **v** will form a parallelogram.
- By transforming it with the rotation matrix $$$A_\theta$$$ we get an analogous set of vectors $$$\\{\mathbf x, \mathbf y\\} \ni \bigg(\mathbf x = \begin{bmatrix} x_1 \\\ x_2 \end{bmatrix}, \mathbf y = \begin{bmatrix} y_1 \\\ 0 \end{bmatrix}\bigg)$$$. The determinant of the matrix of this set of vectors is the area of the parallelogram $$$\mathrm{Area_p} = y_1|x_2|$$$ in $$$\mathbb R^2$$$. Analogously, any n-dimensional parallelogram-like shape in n-space is called a **parallelopiped**. 

We can transform a $$$3\times 3$$$ matrix $$$A$$$ into an upper triangular matrix $$$U$$$. Then $$$\det(A) = \det(U)$$$.

#### Thm 3.3

Let $$$A$$$ be an $$$n\times n$$$ matrix. Then,

1. If $$$B$$$ is obtained from $$$A$$$ via a row swap, $$$\det(B) = \det(A)$$$.


Elementary row operations and det.:

1. scaling: $$$kr_i\to r_i$$$
2. swap: $$$r_1\leftrightarrow r_2$$$
3. $$$kr_i + r_j\to r_j$$$

$$\det B\text{ along row } i = \sum\_{l=1}^n b\_{il}c\_{il} = \sum\_{l=1}^n ka\_{il}c\_{il} = k\det A$$

#### Thm 3.4

Let $$$A,B$$$ be $$$n\times n$$$ matrices. Then,

1. $$$\exists \,A^{-1} \iff \det(A)\neq 0$$$
2. $$$\det(AB) = \det(A)\det(B)$$$
3. $$$\det(A^T)=\det(A)$$$
4. $$$\exists \,A^{-1} \implies \det(A^{-1}) = {1\over \det(A)}$$$

$$$A$$$ is invertible iff $$$A=E_1E_2E_3\dotso E_k$$$
\\[\begin{align}
\det(A) &= \det(E_1\dotso E_k) \\\
&= \det(E_1)\dotso \det(E_k) \\\
&= r \neq 0 \\\
\end{align}\\]

##### Example:

$$A=\begin{bmatrix} 1 & -1 & 2 \\\ -1 & 0 & c \\\ 2 & 1 & 7 \end{bmatrix}$$

Determine values of $$$c$$$ s.t. $$$A^{-1}$$$ does *not* exist

#### Thm 3.5 (Cramer's Rule)

Let $$$A$$$ be an invertible $$$n\times n$$$ matrix, $$$\mathbf b\in \mathbb R^n$$$ and we define $$$M_i = \[\mathbf a_1, \dotsc, \mathbf a\_{i-1}, \mathbf b, \mathbf a\_{i+1},\dotsc,\mathbf a_n\]$$$. Consider $$$A\mathbf x = \mathbf b$$$. Then if $$$x_i={\det(M)_i\over \det(A)}$$$ the vector **x** represents the unique solution.
**Note:** Don't use this on any matrices larger than $$$2\times 2$$$, as it requires $$$n!$$$ operations.

**General Statment:** After a row swap, $$$r_i\leftrightarrow r_j$$$:

- the signs agree iff $$$|i-j|$$$ is even
- $$$B\_{jl}$$$ differs from $$$A\_{il}$$$ by an odd number of swaps of the form $$$r_k\leftrightarrow r\_{k+1} \iff |i-j|$$$ is even

---

### 7/25/14
## $$$\S$$$ 4.1 - [Subspaces](id:anchor6.1)

Suppose $$$\mathbf u, \mathbf v$$$ satisfy the homogeneous equation.
- The set of solutions to $$$A\mathbf x = \mathbf 0$$$ is *closed* under vector addition, scalar multiplication, and contains $$$\mathbf 0$$$.

**Definition:** A set $$$V\subseteq \mathbb R^n$$$ is a **subspace** if it satisfies the following:
1. closure under vector addition
2. closure under scalar multiplication
3. contains $$$\mathbf 0$$$

**Ex:** Show that $$$S=\\{\mathbf x: A\mathbf x = 0\\}$$$

1. Closure under addition: Let $$$\mathbf u,\mathbf v \in S$$$. $$A(\mathbf u + \mathbf v) = A\mathbf u + A\mathbf v = \mathbf 0 + \mathbf 0, \mathbf u + \mathbf v \in S$$
2. Closure under scalar multiples. Let $$$\mathbf u\in S, r\in \mathbb R$$$. $$A(r\mathbf u) = r(A\mathbf u) = r\,\mathbf 0 = \mathbf 0, r\,\mathbf u\in S$$

##### Examples of Subspaces

1. All of $$$\mathbb R^n$$$ (Trivial space).
2. $$$\\{\mathbf 0\\}\subseteq \mathbb R^n$$$ (The zero subspace).
3. Given $$$S=\\{\mathbf v_1, \mathbf v_2,\dotsc,\mathbf v_k\\}\subseteq \mathbb R^n$$$, $$$\mathrm{Span}(S)\subseteq \mathbb R^n$$$ is also a subspace.
	4. Let $$$\mathbf u, \mathbf w\in \mathrm{Span}(S)$$$.
\\[\begin{aligned}
\mathbf u &= \sum\_{i=1}^k c_i\mathbf v_i \\\
\mathbf w &= \sum\_{i=1}^k b_i\mathbf v_i \\\
\mathbf u + \mathbf w &= \sum\_{i=1}^k (b_i + c_i) \mathbf v_i \\\
\therefore \mathbf u + \mathbf v &\subseteq \mathbb R^n \\\
\end{aligned}\\]
	5. closed under scalar multiplication
	6. closed under inclusion of the zero vector.

**Definition:** The **nullspace** of a matrix $$$A$$$ is the set of vectors that satisfy the homogeneous equation.

**Defintion:** The **column space** of a matrix $$$A$$$ is defined as
\\[\begin{align}
\mathrm{Col}(A) &= \\{\mathbf b\in \mathbb R^m: A\mathbf x = \mathbf b\text{ is consistent}\\} \\\
&= \mathrm{span}(\\{\mathbf a_1,\dotsc,\mathbf a_m\\}) \\\
\end{align}\\]
which essentially is the set of all possible linear combinations of the column vectors of $$$A$$$.

## $$$\S$$$ 4.2 - [Bases](id:anchor4.2)

**Definition:** Let $$$V\subseteq \mathbb R^n$$$ be a non-zero subspace. A **basis** is a linearly independent generating set.
- $$$\\{\mathbf e_j\\}$$$ (the standard vectors) form a basis for $$$\mathbb R^n$$$.

#### Thm 4.3

Let $$$S$$$ be a *finite* generating set for some subspace $$$V$$$. There is a subset $$$T\subseteq S$$$ which is a *basis* for $$$V$$$.

For $$$\mathbb R^n$$$:
1. Any generating set contains *at least n* elements.
2. Any linearly indpendent subset contains *at most n* elements.
3. Any basis of $$$\mathbb R^n$$$ has *exactly n* elements.

#### Thm 4.4

Let $$$S$$$ be a linearly independent subset of some subspace $$$V\subseteq \mathbb R^n$$$. Then $$$S$$$ can be extended to a basis of $$$V$$$.
**Proof:**

- $$$S=\\{\mathbf u_1,\dotsc, \mathbf u_k\\}$$$
- $$$V \\ \mathrm{span}(S)\neq \emptyset$$$ provided that $$$S$$$ is not a basis.
- Let $$$\mathbf v\in V \\ \mathrm{span}(S)$$$.
- By previous thm, $$$\\{\mathbf u_1,\dotsc,\mathbf u_k,\mathbf v\\}$$$ is also linearly independent.

#### Thm 4.5

Let $$$V\subseteq \mathbb R^n $$$ be a nonzero subspace. then any two bases have the same number of elements.
**Proof:** 

- Let $$$\\{\mathbf u_1, \dotsc,\mathbf u_k\\}, \\{\mathbf v_1,\mathbf v_p\\}$$$ be two bases for $$$V$$$
- Define $$$A=\[\mathbf u_1, \dotsc, \mathbf u_k\], B=\[\mathbf v_1,\dotsc, \mathbf v_p\]$$$. For each $$$\mathbf v_i, A\mathbf c_i = \mathbf v_i$$$ has a solution. 
- Define $$$C=\[\mathbf c_1,\dotsc,\mathbf c_p\]. C$$$ is a $$$k\times p$$$ matrix.
- $$$AC=B$$$
- $$$B\mathbf x = AC\mathbf x = \mathbf 0$$$
- $$$\forall \mathbf x\in \mathrm{Null}(C) = \\{\mathbf 0\\}$$$
- Therefore, the columns of $$$C$$$ are LI. $$$\therefore p\leq k$$$.
- By a symmetric argument, $$$k\leq p$$$.

**Definition:** The **dimension** of a subspace is the size of one of its bases.

## $$$\S$$$ 4.3 - Dimension

- Col(A) is a subspace. What is its dimension?
	- rank(A).
- Null(A) is a subspace. What is its dimension?
	- nullity(A)

**Definition:** $$$\mathrm{Row}(A) = \\{\mathbf yA:\mathbf y\in \mathbb R^m\\}$$$
- The dimension of the row space is the # of nonzero rows of the rref.
- The basis of the row space will be the nonzero rows of the rref.

$$W\subseteq V\subseteq \mathbb R^n$$
- $$$W$$$ is a subspace of $$$V$$$
- $$$V$$$ is a subspace of $$$\mathbb R^n$$$
- $$$\mathrm{dim}(W)\leq \mathrm{dim}(V)$$$
- If $$$\dim(w)=\dim(V), W=V$$$


--- 

### 7/28/14
## $$$\S$$$ 5.1 - [Eigen-things](id:anchor5)

**Definition:** Let $$$A$$$ be an $$$n\times n$$$ matrix, $$$\lambda$$$ be some *nonzero* scalar, $$$\mathbf v$$$ be some *nonzero* vector in $$$\mathbb R^n$$$. If $$$A\mathbf v = \lambda\mathbf v$$$, we call $$$\lambda$$$ and $$$\mathbf v$$$ an **eigenvalue, eigenvector** pair.

**Ex:**

$$B=\begin{bmatrix}3 & 0 & 0 \\\ 0 & 1 & 2 \\\ 0 & 2 & 1 \\\ \end{bmatrix}$$

Problem: Find a basis for the subspace of vectors with eigenvalue 3.

$$\left\[B - \lambda I: \lambda = 3\right\] = \begin{bmatrix} 0 & 0 & 0 & \\\ 0 & -2 & 2 \\\ 0 & 2 & -2 \\\ \end{bmatrix} = \begin{bmatrix} 0 & 1 & - 1 \\\ 0 & 0 & 0 \\\ 0 & 0 & 0 \\\ \end{bmatrix}$$

A basis for the 3-eigenspace is $$$\left\\{\,\begin{bmatrix} 1 \\\ 0 \\\ 0 \\\ \end{bmatrix}, \begin{bmatrix} 0 \\\ 1 \\\ 1 \\\ \end{bmatrix}\,\right\\}$$$

**Definition:** The subspace of vectors assoc. with the solution $$$(A-\lambda I)\mathbf v = \mathbf 0$$$ is called the **$$$\lambda$$$-eigenspace**.

## $$$\S$$$ 5.2 - [Characteristic Polynomial](id:anchor5.2)

- Consider $$$\det (A-\lambda I)=0$$$ for some $$$n\times n$$$ matrix $$$A$$$ and nonzero scalar $$$\lambda$$$.
	- The eigenvalues of $$$A$$$ are precisely those that satifsy the above equation.
- So what is $$$\det (A-\lambda I)$$$ as a function of $$$\lambda$$$?
	- It's a polynomial of degree $$$n$$$!

**Ex:**

\\[\begin{aligned}
A\, &= \begin{bmatrix}-4 & -3 \\\ 3 & 6 \end{bmatrix} \\\
A - \lambda I \;&= \begin{bmatrix} -4-\lambda & -3 \\\ 3 & 6-\lambda \end{bmatrix} \\\
\det(A-\lambda I)\;&= (4-\lambda)(6-\lambda) - -3(3) \\\
&= -24 -2\lambda + \lambda^2 + 9 \\\
&= \lambda^2 - 2\lambda - 15 = 0 \\\
&= (\lambda - 5)(\lambda + 3) &= 0 \\\
\lambda \;&= -3, 5 \\\
\text{Solve for the general solutions of each }&\lambda\text{-eigenspace to find their respective bases.} \\\
\end{aligned}\\]

	Imaginary eigenvalues tell you the degree of rotation of a matrix
	
**Definition:** The eqn $$$\det(a-\lambda I) = 0$$$ is called the **characteristic equation**. The **multiplicity** (or **algebraic multiplicity**) of an eigenvalue is the corresponding root of the char. poly. If it exists, the $$$\dim(\lambda\text{-space})\leq\mathrm{mult}(\lambda)$$$

---

### 7/30/14
## $$$\S$$$ 5.3 - [Diagonalization of Matrices](id:anchor5.3)

#### Thm 5.2

- An $$$n\times n$$$ matrix $$$A$$$ is **diagonizable** iff there is a basis for $$$\mathbb R^n$$$ consisting entirely of eigenvectors.
- Furthermore, if $$$A=PDP^{-1}$$$, then the columns of $$$P$$$ are the elements of the basis of eigenvectors, and the entries of $$$D$$$ are the corresponding eigenvalues.

**Proof:** Suppose that $$$A$$$ is diagonizable. $$$A=PDP^{-1}$$$ Since $$$P$$$ is invertible, its columns form a basis for $$$\mathbb R^n$$$. Consider $$$AP = PD$$$. 
$$A\mathbf p_j = \lambda_j\mathbf p_j$$ 
Suppose that $$$\\{\mathbf p_1,\dotsc,\mathbf p_n\\}$$$ is a basis of eigenvectors. Let $$$\\{\lambda_1,\dotsc,\lambda_n\\}$$$ be their eigenvalues,, and define $$$P$$$ and $$$D$$$ as before.

#### Thm 5.3

Eigenvectors with distinct eigenvalues are linearly independent.

**Proof:** Suppose $$$\\{\mathbf v_1,\dotsc,\mathbf v_k\\}$$$ eigenvectors, each of which has a distinct eigenvalue $$$\lambda_k$$$. Assume this set is linearly dependent.

##### Test for Diagonalization

A matrix $$$n\times n$$$ $$$A$$$ is diagonizable iff the following holds 

1. For each eigen-$$$\lambda$$$ $$$n = \mathrm{rank}(A, B)$$$

---

### 8/11/14

## $$$\S$$$ 6.1 Dot/Inner Products and their Relations

\\[\begin{align}
\text{dot product:}\;\mathbf u \bullet \mathbf v &= \sum_1^n \mathbf u_i \mathbf v_i \\\
\text{inner product:}\;\<\*,\*>\,: &V\times V\to\mathbb F \\\
&\mathbb R^n \times \mathbb R^n\to \mathbb R \\\
\text{norm:}\;\mid\mid \mathbf u\mid\mid &= \sum_1^n\sqrt{u_i^2} \\\
\end{align}\\]

#### Thm 6.1

For vectors, $$$\mathbf{u\;v\;w}$$$:

1. $$$\mid\mid \mathbf u\mid\mid = \mathbf u \bullet \mathbf u$$$
1. $$$\mathbf u \bullet (\mathbf v + \mathbf w) = \mathbf u \bullet \mathbf v + \mathbf u \bullet \mathbf w)$$$
2. $$$(c\,\mathbf u)\bullet \mathbf v = c(\mathbf u \bullet \mathbf v) = (c\,\mathbf v)\bullet \mathbf u$$$
2. $$$\mid\mid c\,\mathbf u\mid\mid\; = \;\mid c\mid\, \mid\mid \mathbf u \mid\mid $$$

### Extension of the Pythagorean Thm

$$$\mathbf u^2 + \mathbf v^2 = \mathbf w^2$$$

### Cauchy-Schwarz Ineqhuality

$$$\mid\mathbf u \bullet \mathbf v \mid \;\leq \; \mid\mid \mathbf u \mid\mid \; \mid\mid \mathbf v \mid\mid$$$

### Triangle Inequality

$$$\mid\mid \mathbf u + \mathbf v\mid\mid \;\leq\; \mid\mid \mathbf u \mid\mid + \mid\mid \mathbf v \mid\mid$$$

---

## $$$\S$$$ 6.2 - Orthogonality

**Def:** A set of vectors is **orthogonal** iff $$$\mathbf u_i \bullet \mathbf u_j = 0,\, i\neq j$$$

- Orthogonal sets are *linearly independent*
- If you have an orthogonal set, you can determine a linear combination for some vector $$$\mathbf u$$$ with the following formula: $$\mathbf u = \sum_1^n {\mathbf u \bullet \mathbf v_i \over \mid\mid\mathbf v_i\mid\mid^2} \mathbf v_i$$

### Gram-Schmidt Process

Let $$$\\{\mathbf u_1,\dotsc,\mathbf u_k\\}$$$ be a basis for a subspace $$$W$$$.

\\[\begin{align}
\mathbf v_1 &= \mathbf u_1 \\\
\mathbf v_2 &= \mathbf u_2 - {\mathbf u_2 \bullet \mathbf v_1 \over \mid\mid\mathbf v_1\mid\mid^2} \mathbf v_1 \\\
\mathbf v_3 &= \mathbf u_3 - {\mathbf u_3 \bullet \mathbf v_1 \over \mid\mid\mathbf v_1\mid\mid^2} \mathbf v_1 - {\mathbf u_3 \bullet \mathbf v_2 \over \mid\mid\mathbf v_2\mid\mid^2} \mathbf v_2 \\\
&\vdots \\\
\mathbf v_k &= \mathbf u_k - \sum\_{i=1}^{k-1} {\mathbf u_k \bullet \mathbf v_i \over \mid\mid\mathbf v_i\mid\mid^2} \mathbf v_i
\end{align}\\]

---

## $$$\S$$$ 6.3 - Orthogonal Projections

\\[\begin{align}
S^{\perp} &= \\{\mathbf v:\mathbf v\bullet \mathbf s_i = 0, \forall s_i\in S\\} \\\
S^{\perp} &= (\mathrm{span}(S))^{\perp} \\\
\end{align}\\]

#### Thm 6.7

\\[\begin{align}
&\text{Let } W\subseteq \mathbb R^n, \\\
&\\{\mathbf v_1,\dotso,\mathbf v_k\\} \text{ be an orthonormal basis of } W \\\\
&\forall\mathbf u\in \mathbb R^n,\;\exists\! \mathbf w \in W,\, \mathbf z\in W^{\perp},\; \ni \mathbf u = \mathbf w + \mathbf z \\\
&\mathbf w = (\mathbf u\bullet \mathbf v_1)\mathbf v_1 + \dotso + (\mathbf u \bullet \mathbf v_k)\mathbf v_k) \\\
&\mathbf w \text{ is the orthogonal projection of } \mathbf u \text{ onto } W \\\
\end{align}\\]

**Fact:** $$$\dim W + \dim(W^{\perp}) = n$$$, the *ambient space* that you're working in (the $$$n$$$ of your $$$\mathbb R^n$$$.

\\[\begin{align}
\dim(\mathrm{Col}(A)) = \mathrm{rank} \\\
\dim(\mathrm{Null}(A^T)) = m - \mathrm{rank} \\\
\end{align}\\]


---

### 8/8/14

## $$$\S$$$ 6.5 Orthogonal Matrices and Operators

- $$$A\mathbf e_i = \mathbf a_i$$$
- $$$A\mathbf x = \mathbf 0$$$
- $$$A\mathbf v = \lambda v$$$
- $$$\mid\mid A\mathbf v \mid\mid = \mid\mid\mathbf v\mid\mid$$$
	- **Ex:** In $$$\mathbb R^2\; A_\theta = \begin{bmatrix}\cos\theta & \sin\theta \\\ -\sin\theta & \cos\theta \end{bmatrix}$$$

- **Def:** Let $$$A$$$ be an $$$n\times n$$$ matrix, and suppose $$$\mid\mid A\mathbf v \mid\mid = \mid\mid\mathbf v\mid\mid$$$ holds for all vectors $$$\mathbf v$$$. In this case, we say $$$A$$$ is an **orthogonal matrix**.
- Suppose $$$A$$$ is orthogonal:
	1. $$$\mid\mid A\mathbf e_i \mid\mid = \mid\mid \mathbf a_i \mid\mid = \mid\mid \mathbf e_i\mid\mid = 1$$$

	\\[\begin{align}
	\mid\mid \mathbf a_i + \mathbf a_j \mid\mid^2 &= \mid\mid A\mathbf e_i + A \mathbf e_j\mid\mid ^2
	&= \mid\mid A(\mathbf e_i + \mathbf e_j)\mid\mid^2
	&= 2
	&= \mid\mid \mathbf e_i\mid \mid ^2 + \mid\mid \mathbf e_j \mid\mid^2	&= \mid\mid \mathbf a_i\mid \mid ^2 + \mid\mid \mathbf a_j \mid\mid^2
	\end{align}\\]
	
#### Thm 6.9

The following are equivalent for an $$$n\times n$$$ matrix $$$A$$$:

1. $$$A$$$ is orthogonal.
2.	$$$A^TA = I$$$
	$$$(A^{-1} = A^T)$$$
3. $$$(A\mathbf u)\bullet (A\mathbf v) = \mathbf u\bullet \mathbf v = \, <A\mathbf u, A\mathbf v >$$$

#### Thm 6.10

Let $$$P,Q$$$ be two orthogonal matrices

1. $$$PQ$$$ is also orthogonal.
2. $$$Q^{-1}$$$ is also orthogonal.
3. $$$Q^T$$$ ""
4. $$$\det(Q)=\pm 1$$$

**Proof:**
\\[\begin{align}
1 = \det(I) &= \det(QQ^T) \\\
&=\det(Q)\det(Q^T) \\\
&= \det(Q)^2 \\\
\implies &\det(Q) =\pm 1 \\\
\end{align}\\]

**Fact:** 
- Two reflections in a row is a a rotation.
- A reflection followed by a rotation is some reflection.

## $$$\S$$$ 6.6

Let $$$A$$$ be an $$$n\times n$$$ matrix, and suppose that it is diagonalizable. $$$A=PDP^{-1}$$$. If $$$P$$$'s columns are an orthonormal basis, then 
\\[\begin{align}
A=PDP^T \\\
A^T = (PDP^T)^T=P^{T^T}D^TP^T = PDP^T=A \\\
\therefore A\text{ is symmetric.} \\\
\end{align}\\]

#### Thm 6.14

If $$$\mathbf u,\mathbf v$$$ are two eigenvectors of a symmetric matrix $$$A$$$ with different eigenvalues, say $$$\lambda, \mu$$$, then $$$\mathbf u \bullet \mathbf v = 0$$$.
\\[\begin{align}
(A\mathbf u)\bullet \mathbf v &= (\lambda \mathbf u) \mathbf v= \lambda(\mathbf u \bullet \mathbf v) \\\
&= \mathbf u \bullet (A^T\mathbf v) \\\
&= \mathbf u \bullet (A\mathbf v) \\\
&= \mathbf u \bullet (\mu\mathbf v) \\\
&= \mu(\mathbf u \bullet \mathbf v) \\\
\end{align}\\]

#### Thm 6.15

An $$$n\times n$$$ matrix $$$A$$$ is symmetric, iff there is an orthonormal basis of eigenvectors of $$$A$$$ in $$$\mathbb R^n$$$, and in this case, there is an orthogonal matrix and a diagonal matrix $$$\ni A=PDP^T$$$.