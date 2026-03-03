import math
import time
import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from io import BytesIO
from PIL import Image, ImageDraw
from itertools import permutations,combinations
from streamlit_image_coordinates import streamlit_image_coordinates

#pip install xlsxwriter 

#pip install streamlit-image-coordinates
#python -m streamlit run 


def horas_dec_para_relog(horas_dec:float) -> str:#transformando as horas decimais em horas 00:00:00 para os graf interativos

  #tempo é uma lista com seus tempos em horas decimais

    hora = str(math.floor(horas_dec))

    if len(hora) == 1:

        hora = "0" + hora

    minut = (horas_dec - math.floor(horas_dec))*60

    minu = str(math.floor(minut))

    if len(minu) == 1:

        minu = "0" + minu

    segud = (minut - math.floor(minut))*60

    segu = str(int(segud))

    if len(segu) > 2:
        segu = segu[:2]

    if len(segu) == 1:

        segu = "0" + segu

    horas_dec = f"{hora}:{minu}:{segu}"

    return(horas_dec)


def cal_erro(lista1,gx,vx,vy):
    erropor = []

    ytra = np.mean(vy)

    SQresr = 0
    SQexpr = 0

    MAE = 0
    MAPE = 0

    for i in lista1:

        x = vx[i]

        y = eval(gx)

        yr = vy[i]

        if x == 0 or yr == 0:
            pr = 0
        else:
            pr = (abs(yr-y))*100/(yr)

        erropor.append(pr)


        erro_valor = y - yr

        SQresr += (erro_valor) ** 2
        
        SQexpr += (y - ytra) ** 2

        #MAPE e MAE

        if yr != 0:
            erro_porc_v = 100*abs(erro_valor/yr)
        else:
            erro_porc_v = 0
        
        MAE += erro_valor
        MAPE += erro_porc_v



    R2 = SQexpr/(SQexpr+SQresr)

    Erropor = max(erropor)

    return(Erropor,R2,MAE,MAPE)

def ajuste(n:int,
           dict_decisao: dict,
           vx:list,
           vy:list,
           ) -> tuple[str,float,float,float,float,list,list]:
    
    p = len(vx)
    
    gz=[]
    gx = ""
    lista1=np.arange(0,p,1)

    A = np.zeros((n,n))
    B = np.zeros((1,n))[0]
    
    #criação do sistema de equações para econtrar o resultado do ajuste de curva
    for i in lista1:
        
        gs=[]
        x = vx[i]
        
        #enchendo a matriz de x**j
        for j in range(n):
           gs.append(x**j)
        
        #print(gs)
        #trocando os x**j por funções escolhidas
        for j in dict_decisao:

            if dict_decisao[j] != 0:
        
                gs[dict_decisao[j]] = eval(j)
        
        #fazendo as outras partes da matriz
        for k in range(len(A)):

            for j in range(int(len(A))):

                A[j,k] = A[j,k]+(gs[j]*gs[k])
        
        #fazendo a outra matriz, de reusltados
        for j in range(len(A)):
      
            B[j] = B[j] + (vy[i]*gs[j])

    
    #resultado do sistema de quações
    X=np.linalg.inv(A).dot(B)

    
    
    #enchendo a equação com x**j primeiramente
    for j in range(n):

        gz.append("x"+"**"+str(j))
    
    #pegando as funções que foram usadas nessa equação e trocando pelos x**j
    for j in dict_decisao:
        if dict_decisao[j] != 0:
            gz[dict_decisao[j]] = j


    # print(list(X))
    # print(gz)

    #organizando a equação
    for j in range(n):
        
        if X.item(j) < 0:
            
            mais = " "
        
        else:
            
            mais = " +"
        
        if j == 0:
            
            gx = str(gx) +mais+ str(X.item(j))
            
        else:
            gx = str(gx) +mais+ (str(X.item(j))+"*"+str(gz[j]))
        
    # print('Função encontrada:',gx)
    gx=str(gx)
    
    #calculo do erro
    er,R2,MAE,MAPE = cal_erro(lista1,gx,vx,vy)

    return(gx,er,R2,MAE,MAPE,list(X),gz)

def generate_with_unique_nonzero_at_most_k(base_dict, P, k, *, allowed_nonzero=None, values_desc=False):
    """
    Yield dicts with at most k positions non-zero, and all non-zero values are unique.
    """
    keys = list(base_dict.keys())
    n = len(keys)
    pool = list(range(1, P + 1)) if allowed_nonzero is None else list(allowed_nonzero)
    pool.sort(reverse=values_desc)

    # i = number of non-zero positions used (0..k)
    for i in range(0, min(k, n, len(pool)) + 1):
        for nz_positions in combinations(range(n), i):
            for values_perm in permutations(pool, i):
                vals = [0] * n
                for pos, v in zip(nz_positions, values_perm):
                    vals[pos] = v
                yield dict(zip(keys, vals))

def integral_normal(
        min_v:float,
        max_v:float,
        coeficientes:list,
        componentes:list
):  
    eq_final = ""


    for i in range(len(coeficientes)):

        coeficientes[i] = coeficientes[i] / (i+1)

        componentes[i] = componentes[i][:3] + str(int(componentes[i][3:])+1)

        if coeficientes[i] > 0:
            eq_final += f"+{coeficientes[i]}*{componentes[i]}"
        else:
            eq_final += f"{coeficientes[i]}*{componentes[i]}"

    #print(eq_final)

    x = max_v

    valor_cima = eval(eq_final)

    x = min_v

    valor_baixo = eval(eq_final)

    integral = valor_cima - valor_baixo

    return(integral)

def integral_brute_force(
        min_v:list,
        max_v:list,
        equacao:str,
        resolucao:int
):
    
    passo = (max_v - min_v) / resolucao

    integral = 0

    vx = []
    vy = []

    for x in np.arange(min_v,max_v+passo,passo):

        y = eval(equacao)

        vx.append(x)
        vy.append(y)

    for i in range(1,len(vx),1):

        h_tri = vy[i-1] - vy[i]
        
        base  =  vx[i] - vx[i-1]
        
        h_ret = min([(vy[i-1]),vy[i]])
        
        a_tri =  (base * h_tri)/2
        
        a_ret = base * h_ret
        
        integral += a_tri + a_ret

    
    return(integral)


def transformacao_pontos2(
    pontos: list,
    lista_min_max_pontos: list,
    lista_min_max_valores: list,
    y_up_increases=False,  # True if data y increases upward (typical plots)
):
    """
    lista_min_max_valores = [x_min_val, x_max_val, y_min_val, y_max_val]
    lista_min_max_pontos = [(xmin_x, xmin_y), (xmax_x, xmax_y), (ymin_x, ymin_y), (ymax_x, ymax_y)]
    pontos: list of (x_px, y_px)
    """

    # Unpack for clarity
    x_min_val, x_max_val, y_min_val, y_max_val = lista_min_max_valores
    (xmin_x, xmin_y), (xmax_x, xmax_y), (ymin_x, ymin_y), (ymax_x, ymax_y) = lista_min_max_pontos

    # Sanity checks (avoid division by zero)
    if xmax_x == xmin_x:
        raise ValueError("x calibration pixel coordinates have the same x (division by zero).")
    if ymax_y == ymin_y:
        raise ValueError("y calibration pixel coordinates have the same y (division by zero).")

    # Slopes (value per pixel)
    conv_x = (x_max_val - x_min_val) / (xmax_x - xmin_x)

    # For y: image y grows downward. If data y grows upward, slope should be negative w.r.t. pixel y.
    # Use ymin_y and ymax_y in the right order to control sign explicitly:
    if y_up_increases:
        # Increasing data y upward → decreasing pixel y downward
        conv_y = (y_max_val - y_min_val) / (ymin_y - ymax_y)  # note reversed order
        y0 = ymin_y
        val0_y = y_min_val
    else:
        # If your plot's y also increases downward (rare), use normal order
        conv_y = (y_max_val - y_min_val) / (ymax_y - ymin_y)
        y0 = ymin_y
        val0_y = y_min_val

    novos_pontos = []
    for x_px, y_px in pontos:
        # X mapping: anchor at x_min pixel
        x_val = (x_px - xmin_x) * conv_x + x_min_val

        # Y mapping: anchor at y_min pixel (as defined above)
        y_val = (y_px - y0) * conv_y + val0_y

        novos_pontos.append([x_val, y_val])

    #st.markdown(novos_pontos[0])

    vx = []
    vy = []

    for i in novos_pontos:

        vx.append(i[0])
        vy.append(i[1])

    return(vx,vy)


def check_de_valores(dict_decisao,vx) -> dict:
    
    if any(a > 600 for a in vx): #limitação do math.exp(x) que não pode ser maior que +ou- 700
                                 #limitação do math.pi**x que não pode ser maior que +ou- 600

        del dict_decisao['np.exp(x)']
        del dict_decisao['np.exp(-x)']
        del dict_decisao['x*np.exp(-x)']
        del dict_decisao['np.exp(-(x**2))']
        
        del dict_decisao['np.pi**x']
        del dict_decisao['np.sinh(x)']
        del dict_decisao['np.cosh(x)']
        del dict_decisao['np.tanh(x)']

        del dict_decisao['2**x']
        del dict_decisao['3**x']
        del dict_decisao['1.5**x']
        
            
    if min(vx) < 0: #se x < 0, raiz n funciona
        
        del dict_decisao['x**(1/2)']
        del dict_decisao['x**(1/3)']
        del dict_decisao['x**(-1/2)']
            
    if max(vx) > 1.5 : #se em grau > 90, não tem pq usar sen cos e tg

    #     del dict_decisao['np.sin(x)']
    #     del dict_decisao['np.cos(x)']
        del dict_decisao['np.tan(x)']
    #     del dict_decisao['x*np.sin(x)']
    #     del dict_decisao['x*np.cos(x)']
    

    # if any(a == 0 for a in vx): #se o valor de x é muito pequeno, 1/x pode dar valores muitos estranhos

    #     del dict_decisao['x**(-1)']
    #     del dict_decisao['x**(-2)']
    #     del dict_decisao['x**(-3)']
    #     del dict_decisao['x**(-4)']
        
                
    if any(a == 0 for a in vx) or min(vx) < 0: #verifica se algum valor é igual a 0 e se x < 0

        del dict_decisao['np.log(x)']
        del dict_decisao['np.log10(x)']
        del dict_decisao["np.log2(x)"]
    
        
    return(dict_decisao)


def ajuste_de_curva_completo(
        n_max:int,
        todos_dicts:list[dict],
        R2_min:float,
        MAPE_max:float,
        Error_max:float
    ):
    

    lista_de_erros = []
    lista_de_R2 = []
    lista_de_MAPE = []
    lista_de_funcoes = []

    best = []

    texto = st.empty()

    for i in todos_dicts:

        achei = False

        try:
            resultados = ajuste(n_max,i,vx,vy)
            equacao = resultados[0]
            erromax = resultados[1]
            R2      = resultados[2]
            # MAE     = resultados[3]
            MAPE    = resultados[4]
            # coefs   = resultados[5]
            # comps   = resultados[6]

            if isinstance(R2, complex):
                continue
            else:
                lista_de_erros.append(erromax)
                lista_de_R2.append(R2)
                lista_de_MAPE.append(MAPE)
                lista_de_funcoes.append(equacao)

        except:
            pass

        porcen = "Loading... " + str(round(100*(len(lista_de_funcoes)/len(todos_dicts)),2)) + " %"

        texto.markdown(porcen)
        if isinstance(R2, complex) == False:
            if erromax <= Error_max and R2 >= R2_min and MAPE <= MAPE_max:
                achei = True
                st.markdown("Found an equation that fits all criteria:")
                #st.markdown(f"{equacao} {R2 = }, {erromax = }, {MAPE = }")

                best.append(
                    [equacao,R2,erromax,MAPE]
                )

                break
    
    texto.markdown(f"Finished! A total of {len(lista_de_funcoes)} equations were tested.")

    if achei == False:

        st.markdown("Could not find an equation that fits all criteria")

        st.markdown("Here are some of the best results:")

        menor_erro_indx = lista_de_erros.index(min(lista_de_erros))

        best.append(
            [
                lista_de_funcoes[menor_erro_indx],
                lista_de_R2[menor_erro_indx],
                lista_de_erros[menor_erro_indx],
                lista_de_MAPE[menor_erro_indx],
                ]
        )


        menor_MAPE = lista_de_erros.index(min(lista_de_erros))

        best.append(
            [
                lista_de_funcoes[menor_MAPE],
                lista_de_R2[menor_MAPE],
                lista_de_erros[menor_MAPE],
                lista_de_MAPE[menor_MAPE],
                ]
        )

        maior_R2_indx   = lista_de_R2.index(max(lista_de_R2))


        best.append(
            [
                lista_de_funcoes[maior_R2_indx],
                lista_de_R2[maior_R2_indx],
                lista_de_erros[maior_R2_indx],
                lista_de_MAPE[maior_R2_indx],
                ]
        )

    return(best)

def grafico(a_min:float,
            b_max:float,
            resol:int,
            equacao:str,
            nome:str
            )-> go.Figure:

    vxx = []
    vyy = []

    passo = abs(a_min - b_max) / resol

    dx = np.arange(a_min,b_max + passo,passo)

    for x in dx:

        vxx.append(x)
        vyy.append(eval(equacao))

    fig = go.Figure()

    fig.add_trace(go.Scatter(x = vx,y = vy,name = "Points",mode='markers', marker_color='rgba(255, 40, 0, 1)')) # Customize markers

    fig.add_trace(go.Scatter(x = vxx,y = vyy,name=nome,mode='lines',marker_color='rgba(0, 80, 200, 1)')) # Customize line))))

    return(fig)





st.set_page_config("Aplicativo Ajuste de Curva",layout="wide")

st.title("Plot Points with Total Least Squares Method and Integral Calculation")

uploaded_file = st.file_uploader(
    "Image", type=["jpg", "png"],
)


if uploaded_file != None:

    # ---------- Helpers ----------
    # ---- Helpers ----
    def ellipse_coords(p, r=4):
        x, y = p
        return (x-r, y-r, x+r, y+r)

    RED = (255, 0, 0, 255)
    BLUE = (0, 120, 255, 255)

    # ---- Init state ----
    if "points" not in st.session_state:
        st.session_state.points = []             # normal points (red)

    if "calib" not in st.session_state:
        st.session_state.calib = {               # calibration points (blue)
            "x_min": None,
            "x_max": None,
            "y_min": None,
            "y_max": None,
        }

    if "mode" not in st.session_state:
        st.session_state.mode = None             # active calibration target

    if "last_click_sig" not in st.session_state:
        st.session_state.last_click_sig = None   # prevents double event per click
    
    if "last_raw_click" not in st.session_state:
        st.session_state.last_raw_click = None


    # ---- Calibration Buttons ----
    st.subheader("Calibration")

    col = st.columns((2,2,2,2),gap="small")

    x_min = col[0].number_input("X min",format = "%0.5f")
    x_max = col[1].number_input("X max",format = "%0.5f")
    y_min = col[2].number_input("Y min",format = "%0.5f")
    y_max = col[3].number_input("Y max",format = "%0.5f")

    cc1, cc2, cc3, cc4 = st.columns(4)
    if cc1.button("Set x_min"):
        st.session_state.mode = "x_min"
    if cc2.button("Set x_max"):
        st.session_state.mode = "x_max"
    if cc3.button("Set y_min"):
        st.session_state.mode = "y_min"
    if cc4.button("Set y_max"):
        st.session_state.mode = "y_max"

    # ---- Calibration Reset ----
    if st.button("Reset Calibration"):
        st.session_state.calib = {
            "x_min": None,
            "x_max": None,
            "y_min": None,
            "y_max": None,
        }
        st.session_state.mode = None
        st.rerun()


    # ---- Info message ----
    if st.session_state.mode:
        st.info(f"Next click will set **{st.session_state.mode}** (blue).")
    else:
        st.caption("Click to add normal points (red).")


    # ---- Undo / Reset normal points ----
    c1, c2 = st.columns(2)
    with c1:
        if st.button("Undo normal point"):
            if st.session_state.points:
                st.session_state.points.pop()
            st.session_state.last_click_sig = None
            st.rerun()

    with c2:
        if st.button("Reset normal points"):
            st.session_state.points = []
            st.session_state.last_click_sig = None
            st.rerun()


    # ---- Build image ----
    base = Image.open(uploaded_file).convert("RGBA")
    img = base.copy()
    draw = ImageDraw.Draw(img)

    # Draw calibration points in blue
    for key, p in st.session_state.calib.items():
        if p is not None:
            draw.ellipse(ellipse_coords(p), fill=BLUE)

    # Draw normal points in red
    for p in st.session_state.points:
        draw.ellipse(ellipse_coords(p), fill=RED)

    click = None

    # ---- Clickable widget ----
    click = streamlit_image_coordinates(
        img,
        key="clickable_image",
    )


    # ---- Process Click ----

    #print(click)

    if click is not None:

        pt = (click["x"], click["y"])

        # Only process if this is a NEW click position
        if pt != st.session_state.last_raw_click:

            st.session_state.last_raw_click = pt

            if st.session_state.mode is not None:
                # Calibration click
                st.session_state.calib[st.session_state.mode] = pt
                st.session_state.mode = None
            else:
                # Normal point
                st.session_state.points.append(pt)

            st.rerun()
        
        


    # ---- Display calibration values ----
   
    # ---------- Side readout ----------
    with st.expander("Calibration values", expanded=True):
        c1, c2 = st.columns(2)
        with c1:
            
            if st.session_state.calib["x_min"] != None and st.session_state.calib["x_max"] != None:

                st.write("**x_min:**", st.session_state.calib["x_min"][0])
                st.write("**x_max:**", st.session_state.calib["x_max"][0])
        with c2:

            if st.session_state.calib["y_min"] != None and st.session_state.calib["y_max"] != None:

                st.write("**y_min:**", st.session_state.calib["y_min"][1])
                st.write("**y_max:**", st.session_state.calib["y_max"][1])

    st.caption("Legend: blue = calibration points, red = normal points.")


    #st.markdown(st.session_state.points)

    lista_calib = [
        st.session_state.calib["x_min"],
        st.session_state.calib["x_max"],
        st.session_state.calib["y_min"],
        st.session_state.calib["y_max"],
    ]

    lista_calib_valores = [
        x_min,
        x_max,
        y_min,
        y_max
    ]

    #print(st.session_state.points)

    

    
    
    if st.session_state.points != [] and None not in list(dict(st.session_state.calib).values()):

        if x_min != x_max and y_min != y_max:

            vx,vy = transformacao_pontos2(
                st.session_state.points,
                lista_calib,
                lista_calib_valores
                )
            
            dict_data = {
                "X":vx,
                "Y":vy,
            }

            df_points = pd.DataFrame(dict_data)
            
            col_down = st.columns((2,2),gap="medium")
            
            col_down[0].download_button(
                label="Points CSV",
                data=df_points.to_csv(index = False),
                file_name="Points.csv",
                mime="text/csv",
                icon=":material/download:",
            )

            def to_excel_local(df):
                output = BytesIO()
                with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                    df.to_excel(writer, index=False, sheet_name='Sheet1')
                return output.getvalue()

            col_down[1].download_button(
                label="Points Excel(.xlsx)",
                data=to_excel_local(df_points),
                file_name="Points.xlsx",
                mime="xlsx",
                icon=":material/download:",
            )
            
            st.header(r"$\int_{a}^{b}y(x)\, dx$")

            col_ab = st.columns((2,2,2,2),gap="small")
            
            min_local = min(vx)
            max_local = max(vx)

            a_min = col_ab[0].number_input(r"$a$",value=min_local,format="%0.18f")
            b_max = col_ab[1].number_input(r"$b$",value=max_local,format="%0.18f")

            resol = col_ab[2].number_input("Resolution (n of steps)",min_value=1000,max_value=100000,value=1000,step = 1000)

            poly_n = col_ab[3].number_input("Polynomials' Degree",min_value=1,max_value=10,value=1,step=1)

            dict_poly ={
                1:"Linear",
                2:"Quadratic",
                3:"Cubic",
                4:"Quartic",
                5:"Quintic",
                6:"Sextic",
                7:"Septic",
                8:"Octic",
                9:"Notic",
                10:"Decic",
            }


            if len(vx)> 1:
                
                if len(vx) >= poly_n+1:

                    resultado = ajuste(
                        poly_n + 1,
                        {},
                        vx,
                        vy
                    )

                    equacao = resultado[0]
                    erromax = resultado[1]
                    R2   = resultado[2]
                    MAE  = resultado[3]
                    MAPE = resultado[4]
                    coefs = resultado[5]
                    comps = resultado[6]

                    eq_arrumada = ""

                    for i in range(len(coefs)):

                        if coefs[i] < 0:

                            eq_arrumada += fr"{coefs[i]}*{comps[i][0]}^" + "{" +  f"{comps[i][3:]}" + "}"
                        
                        else:
                            eq_arrumada += fr"+{coefs[i]}*{comps[i][0]}^" + "{" +  f"{comps[i][3:]}" + "}"


                    eq_arrumada = eq_arrumada.replace("*x^{1}","*x").replace("*x^{0}","").replace("*",r"\cdot ")
                    
                    integral = integral_normal(
                        a_min,
                        b_max,
                        coefs,
                        comps
                        )

                    "_" * 110

                    st.markdown(fr"$y(x) = {eq_arrumada}$")

                    "_" * 110

                    col_r = st.columns((2,2,2,2,2),gap="small")

                    col_r[0].markdown(r"$\int_{a}^{b}y(x)\, dx =" + f"{integral}$ ")

                    col_r[1].markdown(f"$R^2 = {R2}$")

                    col_r[2].markdown(fr"$MAE = {MAE:.3e}$")

                    col_r[3].markdown(fr"$MAPE = {round(MAPE,2)}\%$")

                    col_r[4].markdown(fr"$Max \,Error = {round(erromax,2)}\%$")

                    "_" * 110

                    fig = grafico(a_min,b_max,resol,equacao,nome = f"{dict_poly[poly_n]} Polynomial")
                
                    st.plotly_chart(fig)


                    

                    dict_decisao = {
                        'x**(2/3)' : 0,
                        'x**(1/2)' : 0,
                        'x**(1/3)' : 0,
                        'x**(-1/2)' : 0,
                        'x**(-1)' : 0,
                        'x**(-2)' : 0,
                        'x**(-3)' : 0,
                        'x**(-4)' : 0,
                        '1.5**x' : 0,
                        '2**x' : 0,
                        '3**x' : 0,
                        'np.pi**x' : 0,
                        'x*(x-1)': 0,
                        'np.sin(x)' : 0,
                        'np.cos(x)' : 0,
                        'np.tan(x)' : 0,
                        #'np.arctan(x)*180/np.pi' : 0,
                        'x*np.sin(x)' : 0,
                        'x*np.cos(x)' : 0,
                        'np.sinh(x)' : 0,
                        'np.cosh(x)' : 0,
                        'np.tanh(x)' : 0,
                        'np.exp(x)' : 0,
                        'np.exp(-x)' : 0,
                        'x*np.exp(-x)' : 0,
                        'np.exp(-(x**2))' : 0,
                        'np.log(x)' : 0,
                        "np.log2(x)" : 0,
                        'np.log10(x)' : 0,
                    }
                    
                    dict_decisao2 = dict(dict_decisao)

                    dict_decisao2 = check_de_valores(dict_decisao2,vx)

                    pode_ir = False

                    with st.form("Parâmetros Ajuste"):

                        col_p = st.columns((2,2,2),gap="small")

                        st.markdown("Input your paramaters to find the best fit curve: (Defaults are based on the curve found above)")

                        poly_n_max = col_p[0].number_input("Polynomials' Degree",min_value=1,max_value=10,value=1,step=1)
                        
                        resol = col_p[1].number_input("Resolution (n of steps)",min_value=1000,max_value=100000,value=10000,step = 1000)


                        substituitions = col_p[2].number_input("Number of Substitutions",
                                                               min_value = 1,
                                                               max_value= poly_n_max+1,
                                                               value = 2,
                                                               help="This value changes the amount of substituitions of $x^n$ for diferent operators, for a value of 2, two values in the equation are going to be switched with values that fit the data. Increasing this number increases computational complexity and lowers speed of the operation.")


                        st.markdown("If you want the best curve possible, leave all values below as 0")

                        

                        R2min = col_p[0].number_input(r"$Minimum \,  R^2$",min_value=0.0,max_value=1.0,value=R2,format="%0.18f")

                        if R2min == 0:
                            R2min = 1

                        MAPEmax = col_p[1].number_input(r"$Maximum \,  MAPE\,[\%]$",value=MAPE,format="%0.18f")

                        erromax = col_p[2].number_input(r"$Maximum \,  Error\,[\%]$",min_value=0.0,value=erromax,format="%0.18f")

                        dict_escolhas = {}

                        col_d = st.columns((2,2,2,2,2,2,2),gap="small")
                        contador = 0

                        for i in dict_decisao:

                            if i not in dict_decisao2:
                                ativado = False
                            else:
                                ativado = True

                            dict_escolhas[i] = col_d[contador].checkbox(i,value=ativado)

                            contador += 1

                            if contador == len(col_d):

                                contador = 0
                        
                        for i in dict_escolhas:

                            if dict_escolhas[i] == False:

                                del dict_decisao[i]
                    

                        #print(dict_escolhas)

                        todos_dicts = [dict_decisao]

                        for i, d in enumerate(generate_with_unique_nonzero_at_most_k(dict_decisao, poly_n_max-1, k=substituitions), start=1):
                            todos_dicts.append(d)

                        qnt_total = len(todos_dicts)

                        sample = 10

                        if qnt_total < 10:
                            sample = int(0.5 * qnt_total)

                        com = time.time()
                        for i in todos_dicts[0:sample]:
                            resultados = ajuste(poly_n_max,i,vx,vy)
                        fim = time.time()

                        tempo_medio = (fim - com) / (10 * 3600)

                        tempo_total = qnt_total * tempo_medio

                        st.markdown(f"Calculation ETA: {horas_dec_para_relog(tempo_total)}")

                        pode_ir = st.form_submit_button('Update values')

                    
                    if st.button("Calculate the best fit curve"):

                        dados = ajuste_de_curva_completo(
                                poly_n_max+1,
                                todos_dicts,
                                R2min,
                                MAPEmax,
                                erromax,
                            )
                        
                        if len(dados) == 1:

                            equacao,R2,erromax,MAPE = dados[0]

                            integral_bruta = integral_brute_force(
                                a_min,
                                b_max,
                                equacao,
                                resol
                            )

                            st.markdown(fr"$y(x) = {equacao}$")

                            col_r2 = st.columns((2,2,2,2,2),gap="small")

                            col_r2[0].markdown(r"$\int_{a}^{b}y(x)\, dx \approx" + f"{integral_bruta}$ ",help = "Since the equation can have any values and any level of complexity, this integral is calculated using the Riemann sum method(using the chosen Resolution), so it is just an approximation")

                            col_r2[1].markdown(f"$R^2 = {R2}$")

                            col_r2[2].markdown(fr"$MAPE = {round(MAPE,2)}\%$")

                            col_r2[3].markdown(fr"$Max \,Error = {round(erromax,2)}\%$")

                            fig2 = grafico(a_min,b_max,resol,equacao,nome = f"Best fit curve")
                
                            st.plotly_chart(fig2,key = 20)
                        
                        else:

                            revolving = [r"lowest $Max \,Error$",r"lowest $MAPE$",r"highest $R^2$"]

                            for i in range(len(dados)):

                                st.markdown(f"Here is the equation with the {revolving[i]}:")

                                equacao,R2,erromax,MAPE = dados[i]

                                integral_bruta = integral_brute_force(
                                    a_min,
                                    b_max,
                                    equacao,
                                    resol
                                )

                                st.markdown(fr"$y(x) = {equacao}$")

                                col_r3 = st.columns((2,2,2,2,2),gap="small")

                                col_r3[0].markdown(r"$\int_{a}^{b}y(x)\, dx \approx" + f"{integral_bruta}$ ",help = "Since the equation can have any values and any level of complexity, this integral is calculated using the Riemann sum method(using the chosen Resolution), so it is just an approximation")

                                col_r3[1].markdown(f"$R^2 = {R2}$")

                                col_r3[2].markdown(fr"$MAPE = {round(MAPE,2)}\%$")

                                col_r3[3].markdown(fr"$Max \,Error = {round(erromax,2)}\%$")

                                fig3 = grafico(a_min,b_max,resol,equacao,nome = f"Best fit curve")
                    
                                st.plotly_chart(fig3,key = i)

                        #print(dados)
                    #print(dict_decisao)

                else:

                    st.markdown(f"The number of points ({len(vx)}) is insufficient for the calculation of the Polynomial of the chosen degree ({poly_n}). Minimum of {poly_n+1} points")

                    st.markdown(f"Try adding more points or lower the Polynomials' Degree")


        else:

            st.markdown(f"The Calibration real values of X min, X max, Y min and Y max have same values that are equal, this should not be possible, please review these values and change them on the top of the page")

            st.markdown(f"{x_min = }")
            st.markdown(f"{x_max = }")
            st.markdown(f"{y_min = }")
            st.markdown(f"{y_max = }")