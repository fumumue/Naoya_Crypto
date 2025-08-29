
unsigned short gcd(unsigned short a, unsigned short b)
{
  int r, tmp;

  if(b==0 || a==0)
  return 0;
  //* 自然数 a > b を確認・入替 
  if (a < b)
  {
    tmp = a;
    a = b;
    b = tmp;
  }

    //b=N;
  //* ユークリッドの互除法 
  r = (a % b);
  while (r != 0)
  {
    a = b;
    b = r;
    r = (a % b);
  }

  //* 最大公約数を出力 
  //printf("最大公約数 = %d\n", b);
  if(b!=1)
  {
    printf("(a,b)!=1\n");
    return -1;
  }

  return b;
}


// 整数の逆数
short inv(short a, short n)
{

    if(gcd(a,n)!=1){
        printf("(a,n)!=1\n");
        return -1;
    }
    short d = n;
    short x = 0;
    short s = 1;
    while (a != 0)
    {
        short q = d / a;
        short r = d % a;
        d = a;
        a = r;
        short t = x - q * s;
        x = s;
        s = t;
    }
    short gcd = d; // $\gcd(a, n)$

    return ((x + n) % (n / d));
}



// 多項式の次数(default)
int deg(vec a)
{
    int i, n = 0, flg = 0;

    // #pragma omp parallel for
    for (i = 0; i < DEG; i++)
    {
        if (a.x[i] != 0)
        {
            n = i;
            flg = 1;
        }
    }
    if (flg == 0)
        return 0;

    return n;
}



vec vadd(vec a, vec b)
{
    int i;
    vec c = {0};

    // printf("deg=%d %d\n",deg(a),deg(b));

    for (i = 0; i < DEG; i++)
        c.x[i] = (a.x[i] + b.x[i]) % N;

    return c;
}



// リーディグタームを抽出(default)
oterm vLT(vec f)
{
    int i;
    oterm t = {0};

    // k = deg (o2v (f));
    for (i = 0; i < DEG; i++)
    {
        // printf("a=%d %d\n",f.t[i].a,f.t[i].n);
        if (f.x[i] > 0)
        {
            t.n = i;
            t.a = f.x[i];
        }
    }

    return t;
}


// aに何をかけたらbになるか
 short
equ( short a,  short b)
{
    // for(short i=0;i<N;i++)
    if (b == 0)
        return 0;
    if (a == 1)
        return b;
    if(a==b)
    return 1;

    return (inv(a, N) * b) % N;
}



// 多項式を単行式で割る
oterm vLTdiv(vec f, oterm t)
{
    oterm tt = {0}, s = {
                        0};

    tt = vLT(f);
    if (tt.n < t.n)
    {
        s.n = 0;
        s.a = 0;
    }
    else if (tt.n == t.n)
    {
        s.n = 0;
        s.a = equ(t.a, tt.a);
    }
    else if (tt.n > t.n)
    {
        s.n = tt.n - t.n;
        s.a = equ(t.a, tt.a);
        // printf("%u\n",s.a);
    }
    else if (t.n == 0 && t.a > 0)
    {
        s.a = (tt.a * inv(t.a, N)) % N;
        s.n = tt.n;
    }

    return s;
}



// 多項式を項ずつ掛ける
vec vterml(vec f, oterm t)
{
    // f = conv(f);
    // ssert(op_verify(f));
    int i;
    vec h = {0};

    // f=conv(f);
    // k = deg (o2v(f));

    for (i = 0; i < DEG; i++)
    {
        // h.t[i].n = f.t[i].n + t.n;
        if (f.x[i] > 0)
            h.x[i + t.n] = (f.x[i] * t.a) % N;
    }

    // h = conv(h);
    //  assert(op_verify(h));
    return h;
}



// 20200816:正規化したいところだがうまく行かない
// 多項式の足し算
vec vsub(vec a, vec b)
{
    vec c = {0};
    // int i, j, k, l = 0;
    vec h = {0}, f2 = {0}, g2 = {0};

    for (int i = 0; i < DEG; i++)
    {
        if (a.x[i] >= b.x[i])
            c.x[i] = (a.x[i] - b.x[i]) % N;
        if (a.x[i] < b.x[i])
            c.x[i] = (N + a.x[i] - b.x[i]) % N;
    }

    return c;
}


int mul = 0, mul2 = 0;
vec vmul(vec a, vec b,int R)
{
    int i, j, k, l;
    vec c = {0};

    k = deg(a);
    l = deg(b);

    if(l+k>N){
        printf("blake %d\n",l+k);
        //exit(1);
    }
    i = 0;
    while (i < k + 1)
    {
        for (j = 0; j < l + 1; j++)
        {
            if (a.x[i] > 0)
                c.x[i + j] = (c.x[i + j] + a.x[i] * b.x[j]) % R;
        }
        i++;
    }

    return c;
}


//整数からベクトル型への変換
vec i2v(unsigned int n)
{
    vec v = {0};
    int i = 0;

    while (n > 0)
    {
        v.x[i++] = n % 2;
        n = (n >> 1);
    }

    return v;
}

//ベクトル型から整数への変換
unsigned v2i(vec v)
{
    unsigned long long int d = 0, i, e = 0;

    for (i = 0; i < 32; i++)
    {
        e = v.x[i];
        d ^= (e << i);
    }

    return d;
}


//モニック多項式にする
vec coeff(vec f,  short d)
{
  int i, j, k;
  vec a, b;

  k = deg((f)) + 1;
  for (i = 0; i < k; i++)
    f.x[i] = (f.x[i]*inv(d,N))%N;

  return f;
}


//多項式の商を取る
vec vdiv(vec f, vec g)
{

  int i = 0, j, n, k;
  vec h = {0}, e = {0}, tt = {0};
  oterm a, b = {0}, c = {0};

  if (vLT(f).n == 0 && vLT(g).a == 0)
  {
    printf("baka^\n");
    //return f;
    exit(1);
  }
  if (vLT(g).a == 0)
  {
    exit(1);
  }
  if (vLT(g).n == 0 && vLT(g).a > 1)
    return coeff(f, vLT(g).a);

  k = deg(g);
  b = vLT(g);
  if (b.a == 1 && b.n == 0)
    return f;
  if (b.a == 0 && b.n == 0)
  {
    printf("baka in vdiv\n");
    exit(1);
  }
  if (deg((f)) < deg((g)))
  {
    return f;
    //  a=LT(f);
  }

  i = 0;
  while (vLT(f).n > 0 && vLT(g).n > 0)
  {
    c = vLTdiv(f, b);
    assert(c.n < DEG);
    tt.x[c.n] = c.a;
    //i++;

    h = vterml(g, c);

    f = vsub(f, h);
    if (deg((f)) == 0 || deg((g)) == 0)
    {
      //printf ("blake2\n");
      break;
    }

    if (c.n == 0)
      break;
  }

  // tt は逆順に入ってるので入れ替える
  return tt;
}


int wt(vec e){
    int i,count=0;

    for(i=0;i<N;i++)
    if(e.x[i]>0)
    count++;

    return count;
}


