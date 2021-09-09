
#pasos:
#obtener h_orth: h_orth = AltOrth(p)
#obtener h_geop: h_geop = AltGeop(p) #esta no depende de latlon
#obtener g0 = g(lat, h_orth=0) = Geographic(LatLon(p), AltOrth(0)) |> gravity

#las cosas que pueden cambiar en un ISA model son T0, p0 y la velocidad del
#viento
#pero las ecuaciones para el calculo de p y t son siempre las mismas
#entonces basicamente lo que define el ISA model son T0, p0, wind y
#turbulence, junto con los sistemas dinamicos que definen la variacion
#temporal de esas variables. pero la estructura es la que es. entonces,
#le pedimos al modelo que nos devuelva T0, p0 en una vw y turbulencia en una
#Abstract3DPosition. con T0, p0

#podemos permitir que varien T0 y p0 con el tiempo.

#meteo model tendra:
#ISA model, static o dynamic
#wind model, static o dynamic

#un ISA model proporcionara un method que dara T, p, rho, a, dados loc2d h_geop y t.
#internamente, los calculara con sus valores almacenados de T0 y p0, que podran
#variar o no dependiendo de si es un ISA static o un ISA dynamic. si es un ISA
#dynamic, tendremos que hacerle step, claro

#un wind model proporcionara un method que dara v_ew_n, dados pos3d y t. aqui
#nuevamente tendremos diferentes modelos posibles. un constant se limitara a
#devolver un valor pre-almacenado que podemos tocar haciendo y=u. lo devolvemos
#en una y immutable, para que el deepcopy no nos joda. el siguiente paso de
#complejidad es un static, que basicamente implementa una variacion espacial del
#viento, es decir, un method get_wind_velocity(p::Abstract3DPosition). esto
#puede ser tan sencillo como una variacion parametrica en funcion de la altitud.

#lo que cambia entre unos modelos y otros es el valor de T0 y p0, que es
#funcion de la Abstract2DLocation y de t. esto nos lo da un ISA model que puede
#ser static o dynamic, y tendra siempre un method que, dados loc y t, nos
#devuelva T0 y p0 locales. con eso podemos entrar en las ecs de la ISA y sacar
#p,T,rho,a