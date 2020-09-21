# data format: width, height, pos_x, pos_y
# where pos_x, pos_y marks the bottom-left corner
data = """
.3 .2 .0 .8
.7 .2 .3 .8
.2 .7 .0 .1
.25 .5 .2 .3
.25 .5 .45 .3
.3 .2 .7 .6
.2 .2 .2 .1
.4 .1 .0 .0
.15 .15 .4 .15
.15 .15 .55 .15
.15 .15 .4 .0
.15 .15 .55 .0
.3 .6 .7 .0
"""

class Patch:
    def __init__(self, w, h, x, y):
        self.w = w
        self.h = h
        self.x = x
        self.y = y

data = [Patch(*map(float, x.split())) for x in data.split('\n') if x]

for i, p in enumerate(data):
    with open('t{:03d}.amr'.format(i), 'w') as f:
        f.write("{} {} {} {}\n".format(p.x, p.x + p.w, p.y, p.y + p.h))
        f.write("1\n")
        f.write("0.01 0.01 1\n")
        f.write('1 1 {} {}'.format(int(p.w / 0.01), int(p.h / 0.01)))
