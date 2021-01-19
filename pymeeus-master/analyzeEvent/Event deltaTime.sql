/*
deltaTime Events PyMeeus vs. Thuillot
*/
WITH

/* PyMeeus Events mit Compound-Type als eine Tabelle definieren */
PyMeeus AS
(
SELECT DateTimeEventTT, Satellite || Type1 || Type2 || Type3 AS Type, DateTimeEventTT, 
Satellite, Type1, Type2, Type3
FROM Event
WHERE Origin="PyMeeus"
),

/* Thuillot Events mit Compound-Type als eine weitere Tabelle definieren */
Thuillot AS
(
SELECT DateTimeEventTT, Satellite || Type1 || Type2 || Type3 AS Type, DateTimeEventTT
FROM Event
WHERE Origin="Thuillot"
)

/* pro PyMeeus Event aus jeweils allen identischen Thuillot Event-Typen den mit der kleinsten absoluten zeitlichen
   Abweichung auswählen und pro PyMeeus Event zeitliche Abweichung zum korrespondieren Thuillot Event anzeigen */
SELECT p.DateTimeEventTT, p.Satellite, p.Type1, p.Type2, p.Type3, 
strftime('%s', p.DateTimeEventTT) - strftime('%s', t.DateTimeEventTT) AS deltaTime, 
MIN( ABS( strftime('%s', p.DateTimeEventTT) - strftime('%s', t.DateTimeEventTT) )) AS minDeltaTime
FROM (PyMeeus p INNER JOIN Thuillot t
ON p.Type = t.Type)
GROUP BY p.DateTimeEventTT
HAVING minDeltaTime <= 3600
/* falls minDeltaTime > 1 Stunde : interpretiere als kein passender Event in Thuillot vorhanden */

