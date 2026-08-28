NEWS:
Elkészült a 2 site no external field coupling ág megoldása, a progress_tracker mappában van a latex dokumentum, a code mappában pedig a python script, ami elvégezte a klasszifikációt.
Lefuttattam a kódot többféle megkötés mellett is, a legáltalánosabb ansatz, ami mellett még lefutott, az volt, ahol csak a 00, 0-, -0, +0 tagok voltak megtiltva (érdekes, hogy +0 tiltással 52.942 s, míg 0+ tiltással 45.758 s volt a futási idő, holott ha jól gondolom, ugyanazt a láncot adják a periodikus határfeltétel miatt). Abban az esetben, ahol 00, 0-, -0 tagok voltak csak megtiltva, több óra alatt sem futott le.

