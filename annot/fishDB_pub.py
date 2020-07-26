from flask import Flask, render_template, redirect, url_for, g, request
import pandas as pd
import sqlite3
import datetime
import random

app = Flask(__name__)

#### DB for pattern annotation
# DATABASE = 'fishDB_pub.db'
DATABASE = 'fishDB_pub_empty.db'
TABLE = 'fishdb'

#### Annotator's name
# ANNOTATOR = 'seita'
ANNOTATOR = 'anonymous'

def get_db():
    db = getattr(g, '_database', None)
    if db is None:
        db = g._database = sqlite3.connect(DATABASE)
        db.row_factory = sqlite3.Row

    return db

def query_db(query, args=(), one=False):
    cur = get_db().execute(query, args)
    rv = cur.fetchall()
    cur.close()
    return (rv[0] if rv else None) if one else rv

@app.teardown_appcontext
def close_connection(exception):
    db = getattr(g, '_database', None)
    if db is not None:
        db.close()

#### Dictionary for Japanese name
with app.app_context():
    sql_str = 'SELECT DISTINCT family, jpn_fam FROM fishdb'
    jpn_fams = dict(query_db(sql_str))

    sql_str = 'SELECT DISTINCT species, jpn_name FROM fishdb'
    jpn_names = dict(query_db(sql_str))

#### for time recording
start_time = datetime.datetime.now()
current_sp = 0
sp100_mode = False
SP100_NUM = 100

#### for auto mode
random_fam = False
current_fam = ""

#### where to go back?
next_url = ""

@app.route('/')
def hello():
    return redirect(url_for('show_fam_list'))

@app.route('/fam_list')
def show_fam_list():
    sort = request.args.get('sort', default=0, type=int)
    sql_str = 'SELECT family, jpn_fam, genus, ' \
              'COUNT(DISTINCT (CASE WHEN is_checked=1 THEN species ELSE NULL END)) AS sp_chk_num, ' \
              'COUNT(DISTINCT species) AS sp_num ' \
              'FROM fishdb GROUP BY genus'
    sql_str = 'SELECT family, jpn_fam, ' \
              'COUNT(CASE WHEN sp_chk_num = sp_num THEN genus ELSE NULL END) AS gen_comp_num, ' \
              'COUNT(genus) AS gen_num, ' \
              'SUM(sp_chk_num) AS sp_chk_num, ' \
              'SUM(sp_num) AS sp_num ' \
              'FROM (' + sql_str + ') GROUP BY family'
                
    if (sort==0): ### A-Z
        sql_str += ' ORDER BY family ASC'
    elif (sort==1): ### gen_num ASC
        sql_str += ' ORDER BY gen_num ASC'
    elif (sort==2): ### gen_num DESC
        sql_str += ' ORDER BY gen_num DESC'
    fam_gen = query_db(sql_str)

    sum_gen_comp_num = sum([row['gen_comp_num'] for row in fam_gen])
    sum_gen_num = sum([row['gen_num'] for row in fam_gen])
    sum_sp_chk_num = sum([row['sp_chk_num'] for row in fam_gen])
    sum_sp_num = sum([row['sp_num'] for row in fam_gen])

    return render_template('fam_list.html', 
                           sum_gen_comp_num = sum_gen_comp_num,
                           sum_gen_num = sum_gen_num,
                           sum_sp_chk_num = sum_sp_chk_num,
                           sum_sp_num = sum_sp_num,
                           fam_gen = fam_gen, 
                           sp100_mode = sp100_mode,
                           current_sp = current_sp)

@app.route('/gen_list/<fam>')
def show_gen_list(fam):
    sort = request.args.get('sort', default=0, type=int)
    
    jpn_fam = jpn_fams[fam]

    ### Sometimes the same genus can be found in different families...
    ### Counting all of them based on the genus list
    
    ### get genus list
    sql_str = 'SELECT genus FROM fishdb WHERE family = ? GROUP BY genus'
    gen_list = [row[0] for row in query_db(sql_str, [fam])]

    ### get sp_num
    format_str = ','.join(list('?' * len(gen_list)))
    sql_str = 'SELECT genus, ' \
              'COUNT(DISTINCT species) AS sp_num, ' \
              'COUNT(DISTINCT CASE WHEN is_checked=1 THEN species ELSE NULL END) AS sp_chk_num ' \
              'FROM fishdb WHERE genus IN({}) GROUP BY genus'.format(format_str)

    if (sort==0): ### A-Z
        sql_str += ' ORDER BY genus ASC'
    elif (sort==1): ### gen_num ASC
        sql_str += ' ORDER BY sp_num ASC'
    elif (sort==2): ### gen_num DESC
        sql_str += ' ORDER BY sp_num DESC'

    gen_sp = query_db(sql_str, gen_list)
    return render_template('gen_list.html', fam=fam, jpn_fam=jpn_fam, gen_sp=gen_sp, sp100_mode=sp100_mode, current_sp=current_sp)

@app.route('/sp_list/<gen>')
def show_sp_list(gen):
    ### Sometimes the same species can be found in different families...
    ### Counting all of them
    sql_str = 'SELECT DISTINCT family, jpn_fam FROM fishdb WHERE genus = ? GROUP BY family, jpn_fam'
    fams = query_db(sql_str, [gen])
    
    sql_str = 'SELECT species, jpn_name, ' \
              'COUNT(img_file) AS pic_num, ' \
              'COUNT(CASE WHEN is_checked=1 THEN img_file ELSE NULL END) AS pic_chk_num, ' \
              'COUNT(mono=1 OR NULL) AS mono_num, ' \
              'COUNT(area_fill=1 OR NULL) AS area_fill_num, ' \
              'COUNT(stripe_vert=1 OR NULL) AS stripe_vert_num, ' \
              'COUNT(stripe_horz=1 OR NULL) AS stripe_horz_num, ' \
              'COUNT(stripe_diag=1 OR NULL) AS stripe_diag_num, ' \
              'COUNT(stripe_maze=1 OR NULL) AS stripe_maze_num, ' \
              'COUNT(spot_dark=1 OR NULL) AS spot_dark_num, ' \
              'COUNT(spot_light=1 OR NULL) AS spot_light_num, ' \
              'COUNT(eyespot=1 OR NULL) AS eyespot_num, ' \
              'COUNT(saddle=1 OR NULL) AS saddle_num, ' \
              'COUNT(blotch=1 OR NULL) AS blotch_num, ' \
              'COUNT(no_use=1 OR NULL) AS no_use_num ' \
              'FROM fishdb WHERE genus = ? GROUP BY species, jpn_name'
    sp_pic = query_db(sql_str, [gen])

    return render_template('sp_list.html', fams=fams, gen=gen, sp_pic=sp_pic, sp100_mode=sp100_mode, current_sp=current_sp)

@app.route('/pic_list/<sp>')
def show_pic_list(sp):
    global next_url
    next_url = request.args.get('next_url', default="", type=str) ### where to go back?

    jpn_name = jpn_names[sp]

    gen = sp.split()[0]
    ### Sometimes the same species can be found in different families...
    ### Counting all of them
    sql_str = 'SELECT DISTINCT family, jpn_fam FROM fishdb WHERE genus = ? GROUP BY family, jpn_fam'
    fams = query_db(sql_str, [gen])

    sql_str = 'SELECT img_file, img_file_path, is_checked,' \
              'mono, area_fill, ' \
              'stripe_vert, stripe_horz, stripe_diag, stripe_maze, ' \
              'spot_dark, spot_light, eyespot, saddle, blotch, no_use, ' \
              'datetime(timestamp) as timestamp, annotator ' \
              'FROM fishdb WHERE species = ? ORDER BY is_checked DESC'
    pic_path = query_db(sql_str, [sp])

    if (len(pic_path)==1):  ### only one image ->  directly to annotate_pic
        if next_url == "":
            next_url = url_for('show_sp_list', gen=gen)
        return redirect(url_for('annotate_pic', img_file=pic_path[0]['img_file']))
    else:
        if next_url == "":
            next_url = url_for('show_sp_list', gen=gen)
        return render_template('pic_list.html', fams=fams, gen=gen, sp=sp, jpn_name=jpn_name, pic_path=pic_path, next_url=next_url, sp100_mode=sp100_mode, current_sp=current_sp)

def get_annotation(img_file):
    annotation = get_db().execute(
        'SELECT family, genus, species, jpn_name, db, img_file_path, '
        ' mono, area_fill, '
        ' stripe_vert, stripe_horz, stripe_diag, stripe_maze, '
        ' spot_dark, spot_light, eyespot, saddle, blotch, no_use'
        ' FROM fishdb'
        ' WHERE img_file = ?',
        (img_file,)
    ).fetchone()

    if annotation is None:
        abort(404, "Annotation for file {0} doesn't exist.".format(img_file))

    return annotation

@app.route('/annotate_pic/<img_file>', methods=('GET', 'POST'))
def annotate_pic(img_file):
    global current_sp
    
    print(next_url)

    annotation = get_annotation(img_file)
    
    if request.method == 'POST':
        mono, area_fill = 0, 0
        stripe_vert, stripe_horz, stripe_diag, stripe_maze = 0, 0, 0, 0
        spot_dark, spot_light, eyespot, saddle, blotch, no_use = 0, 0, 0, 0, 0, 0
        is_checked = 0
        
        if request.form.get('mono'):
            mono = 1
        if request.form.get('area_fill'):
            area_fill = 1
        if request.form.get('stripe_vert'):
            stripe_vert = 1
        if request.form.get('stripe_horz'):
            stripe_horz = 1
        if request.form.get('stripe_diag'):
            stripe_diag = 1
        if request.form.get('stripe_maze'):
            stripe_maze = 1
        if request.form.get('spot_dark'):
            spot_dark = 1
        if request.form.get('spot_light'):
            spot_light = 1
        if request.form.get('eyespot'):
            eyespot = 1
        if request.form.get('saddle'):
            saddle = 1
        if request.form.get('blotch'):
            blotch = 1
        if request.form.get('no_use'):
            no_use = 1
        
        if (mono or area_fill or \
            stripe_vert or stripe_horz or stripe_diag or stripe_maze or \
            spot_dark or spot_light or eyespot or saddle or blotch or no_use):
            is_checked = 1
            
        db = get_db()
        now = datetime.datetime.now()
        db.execute(
            'UPDATE fishdb SET '
            ' mono = ?, '
            ' area_fill = ?, '
            ' stripe_vert = ?, '
            ' stripe_horz = ?, '
            ' stripe_diag = ?, '
            ' stripe_maze = ?, '
            ' spot_dark = ?, '
            ' spot_light = ?, '
            ' eyespot = ?, '
            ' saddle = ?, '
            ' blotch = ?, '
            ' no_use = ?, '
            ' is_checked = ?, '
            ' timestamp = ?, '
            ' annotator = ? '
            ' WHERE img_file = ?',
            (mono, area_fill, stripe_vert, stripe_horz, stripe_diag, stripe_maze, \
             spot_dark, spot_light, eyespot, saddle, blotch, no_use, is_checked, \
             now, ANNOTATOR, img_file)
        )
        db.commit()
        
        #### for time recording
        current_sp += 1
        print("current_sp: ", current_sp)
        
        if (sp100_mode and current_sp >= SP100_NUM):
            return redirect(url_for('sp100', command='lap',
                                     next_url=next_url,
                                     fam=[annotation['family']],
                                     gen=[annotation['genus']]))
        
        if "search" in next_url:
            return redirect(next_url+"#"+img_file)
        else:
            return redirect(next_url+"#"+annotation['species'])
    
    else:
        sql_str = 'SELECT DISTINCT family, jpn_fam FROM fishdb WHERE genus = ? GROUP BY family, jpn_fam'
        fams = query_db(sql_str, [annotation['genus']])

        return render_template('annotate_pic.html', fams=fams, img_file=img_file, annotation=annotation, sp100_mode=sp100_mode, current_sp=current_sp)

@app.route('/query')
def query():
    return render_template('query.html')
    
@app.route('/search')
def search():
    global next_url
    
    is_checked = request.args.get('is_checked', default=1, type=int)

    mono = request.args.get('mono', default=-1, type=int)
    area_fill = request.args.get('area_fill', default=-1, type=int)

    stripe_any = request.args.get('stripe_any', default=-1, type=int)
    stripe_dir = request.args.get('stripe_dir', default=-1, type=int)
    stripe_vert = request.args.get('stripe_vert', default=-1, type=int)
    stripe_horz = request.args.get('stripe_horz', default=-1, type=int)
    stripe_diag = request.args.get('stripe_diag', default=-1, type=int)
    stripe_maze = request.args.get('stripe_maze', default=-1, type=int)

    spot_any = request.args.get('spot_any', default=-1, type=int)
    spot_noeye = request.args.get('spot_noeye', default=-1, type=int)
    spot_dark = request.args.get('spot_dark', default=-1, type=int)
    spot_light = request.args.get('spot_light', default=-1, type=int)
    eyespot = request.args.get('eyespot', default=-1, type=int)

    saddle = request.args.get('saddle', default=-1, type=int)
    blotch = request.args.get('blotch', default=-1, type=int)
    no_use = request.args.get('no_use', default=-1, type=int)

    gen = request.args.get('gen', default="", type=str)
    fam = request.args.get('fam', default="", type=str)

    timestamp = request.args.get('timestamp', default="", type=str)
    annotator = request.args.get('annotator', default="", type=str)

    sql_str = 'SELECT family, jpn_fam, genus, species, jpn_name, img_file, img_file_path, ' \
              'mono, area_fill, ' \
              'stripe_vert, stripe_horz, stripe_diag, stripe_maze, ' \
              'spot_dark, spot_light, eyespot, saddle, blotch, no_use ' \
              'FROM fishdb '

    where_str = 'WHERE is_checked=? '
    if (mono == 0):        where_str += 'AND mono = 0 '
    if (mono == 1):        where_str += 'AND mono = 1 '

    if (area_fill == 0):   where_str += 'AND area_fill = 0 '
    if (area_fill == 1):   where_str += 'AND area_fill = 1 '

    if (stripe_any == 0): where_str += 'AND (stripe_vert=0 AND stripe_horz=0 AND stripe_diag=0 AND stripe_maze=0) '
    if (stripe_any == 1): where_str += 'AND (stripe_vert=1 OR stripe_horz=1 OR stripe_diag=1 OR stripe_maze=1) '

    if (stripe_dir == 0): where_str += 'AND (stripe_vert=0 AND stripe_horz=0 AND stripe_diag=0) '
    if (stripe_dir == 1): where_str += 'AND (stripe_vert=1 OR stripe_horz=1 OR stripe_diag=1) '

    if (stripe_vert == 0): where_str += 'AND stripe_vert = 0 '
    if (stripe_vert == 1): where_str += 'AND stripe_vert = 1 '

    if (stripe_horz == 0): where_str += 'AND stripe_horz = 0 '
    if (stripe_horz == 1): where_str += 'AND stripe_horz = 1 '

    if (stripe_diag == 0): where_str += 'AND stripe_diag = 0 '
    if (stripe_diag == 1): where_str += 'AND stripe_diag = 1 '

    if (stripe_maze == 0): where_str += 'AND stripe_maze = 0 '
    if (stripe_maze == 1): where_str += 'AND stripe_maze = 1 '

    if (spot_any == 0):    where_str += 'AND (spot_dark=0 AND spot_light=0 AND eyespot=0 AND saddle=0) '
    if (spot_any == 1):    where_str += 'AND (spot_dark=1 OR spot_light=1 OR eyespot=1 OR saddle=1) '

    if (spot_noeye == 0):  where_str += 'AND (spot_dark=0 AND spot_light=0) '
    if (spot_noeye == 1):  where_str += 'AND (spot_dark=1 OR spot_light=1) '

    if (spot_dark == 0):   where_str += 'AND spot_dark = 0 '
    if (spot_dark == 1):   where_str += 'AND spot_dark = 1 '

    if (spot_light == 0):  where_str += 'AND spot_light = 0 '
    if (spot_light == 1):  where_str += 'AND spot_light = 1 '

    if (eyespot == 0):     where_str += 'AND eyespot = 0 '
    if (eyespot == 1):     where_str += 'AND eyespot = 1 '

    if (saddle == 0):      where_str += 'AND saddle = 0 '
    if (saddle == 1):      where_str += 'AND saddle = 1 '

    if (blotch == 0):      where_str += 'AND blotch = 0 '
    if (blotch == 1):      where_str += 'AND blotch = 1 '

    if (no_use == 0):      where_str += 'AND no_use = 0 '
    if (no_use == 1):      where_str += 'AND no_use = 1 '

    if (fam != ""):        where_str += 'AND family = "' + fam + '" '
    if (gen != ""):        where_str += 'AND genus = "' + gen + '" '
    
    order_str = 'ORDER BY family, genus, species'
    
    search_result = query_db(sql_str + where_str + order_str,
                             [is_checked])

    hit_fam_num = len(set([row['family'] for row in search_result]))
    hit_gen_num = len(set([row['genus'] for row in search_result]))
    hit_sp_num = len(set([row['species'] for row in search_result]))
    
    next_url = request.full_path
    print(next_url)
    
    return render_template('search_results.html', hit_fam_num=hit_fam_num, hit_gen_num=hit_gen_num, hit_sp_num=hit_sp_num,
                           hit_num=len(search_result), search_result=search_result, sp100_mode=sp100_mode, current_sp=current_sp)

@app.route('/status')
def status():
    pass

"""
@app.route('/auto_fam')
def auto_fam():
    global random_fam
    global current_fam

    fam = request.args.get('fam', default="", type=str)
    if (fam == ""):
        random_fam = True

        ### family, gen_comp_num, gen_num (gen_comp_num < gen_num)
        sql_str = 'SELECT family, genus, ' \
                  'COUNT(DISTINCT (CASE WHEN is_checked=1 THEN species ELSE NULL END)) AS sp_chk_num, ' \
                  'COUNT(DISTINCT species) AS sp_num ' \
                  'FROM fishdb GROUP BY genus'
        sql_fam = 'SELECT family, ' \
                  'COUNT(CASE WHEN sp_chk_num = sp_num THEN genus ELSE NULL END) AS gen_comp_num, ' \
                  'COUNT(genus) AS gen_num ' \
                  'FROM (' + sql_str + ') ' \
                  'GROUP BY family ' \
                  'HAVING (gen_comp_num < gen_num)'
        
        fam_gen_num = query_db(sql_fam)

        ### min of gen_comp_num
        gen_comp_num_min = min([row['gen_comp_num'] for row in fam_gen_num])
        print("gen_comp_num_min = ", gen_comp_num_min)
        
        ### list of families to be checked
        fam_to_be_chk = []
        for row in fam_gen_num:
            if ((row['gen_comp_num'] == gen_comp_num_min) and (row['gen_comp_num'] < row['gen_num'])):
                fam_to_be_chk.append(row['family'])

        random_fam = random.choice(fam_to_be_chk)
        print("random_fam = ", random_fam)

        current_fam = random_fam

    else:
        random_fam = False
        current_fam = fam

    return redirect(url_for('auto_gen', fam=current_fam))

    
@app.route('/auto_gen/<fam>')
def auto_gen(fam):
          
    ### get the list of genera within family
    sql_gen = 'SELECT DISTINCT genus FROM fishdb ' \
              'WHERE family="' + fam + '"'
    
    ### genus WHERE (sp_chk_num < sp_num) ORDER BY sp_num
    sql_gen = 'SELECT genus, '\
              'COUNT(DISTINCT (CASE WHEN is_checked=1 THEN species ELSE NULL END)) AS sp_chk_num, ' \
              'COUNT(DISTINCT species) AS sp_num ' \
              'FROM fishdb ' \
              'WHERE (genus in (' + sql_gen + ')) ' \
              'GROUP BY genus ' \
              'HAVING ( sp_chk_num < sp_num ) ' \
              'ORDER BY sp_num DESC'


    # ### family, genus WHERE (sp_chk_num < sp_num) ORDER BY sp_num
    # sql_str = 'SELECT family, genus, ' \
    #           'COUNT(DISTINCT (CASE WHEN is_checked=1 THEN species ELSE NULL END)) AS sp_chk_num, ' \
    #           'COUNT(DISTINCT species) AS sp_num ' \
    #           'FROM fishdb WHERE family="' + random_fam + '" GROUP BY genus'
        
    # sql_gen = 'SELECT genus, sp_chk_num, sp_num ' \
    #           'FROM (' + sql_str + ') ' \
    #           'WHERE (sp_chk_num < sp_num) ' \
    #           'ORDER BY family, sp_num DESC'

    print(sql_gen)

    gen_sp_num = query_db(sql_gen)

    gen_to_be_chk = [row['genus'] for row in gen_sp_num]
    gen_sp_chk_num = [row['sp_chk_num'] for row in gen_sp_num]
    gen_sp_num = [row['sp_num'] for row in gen_sp_num]
    print("gen_to_be_chk:")
    print(gen_to_be_chk)
    print(gen_sp_chk_num)
    print(gen_sp_num)

    if (len(gen_to_be_chk) > 0):
        return redirect(url_for('auto_sp', gen=gen_to_be_chk[0]))
    else:
        if (random_fam):
            return redirect(url_for('auto_fam', fam=""))
        else:
            return redirect(url_for('show_gen_list', fam=fam))

@app.route('/auto_sp/<gen>')
def auto_sp(gen):
    global next_url
    
    sql_str = 'SELECT species, ' \
              'COUNT(CASE WHEN is_checked=1 THEN img_file ELSE NULL END) AS pic_chk_num ' \
              'FROM fishdb WHERE genus = ? GROUP BY species'
    sql_sp = 'SELECT species ' \
             'FROM (' + sql_str + ') '\
             'WHERE (pic_chk_num = 0)'

    sp_not_chk = query_db(sql_sp, [gen])

    sp_to_be_chk = [row['species'] for row in sp_not_chk]

    print(sp_to_be_chk)

    next_url = request.full_path
    if (len(sp_to_be_chk) > 0):
        return redirect(url_for('show_pic_list', sp=sp_to_be_chk[0], next_url=next_url))
    else:
        if (random_fam):
            return redirect(url_for('auto_fam', fam=""))
        else:
            return redirect(url_for('auto_gen', fam=current_fam))


@app.route('/sp100/<command>')
def sp100(command):
    global start_time
    global current_sp
    global sp100_mode

    global next_url
    
    fam = request.args.get('fam', default="", type=str)
    gen = request.args.get('gen', default="", type=str)
    next_url = request.args.get('next_url', default="", type=str)

    if (command == "start"):
        start_time = datetime.datetime.now()
        current_sp = 0
        sp100_mode = True

        print(sp100_mode)
        print(request.referrer)
        
        return redirect(request.referrer)

    if (command == "lap"):
        
        now = datetime.datetime.now()
        record_time = now - start_time
        start_time = datetime.datetime.now()
        current_sp = 0
        
        sql_str = 'INSERT INTO sp100(timestamp, record) VALUES(?, ?)'
        query_db(sql_str, (str(now), str(record_time)))
        get_db().commit()
        
        sql_str = 'SELECT timestamp, record FROM sp100 ORDER BY record'
        sp100 = query_db(sql_str)
        df_sp100 = pd.read_sql(sql_str, get_db())

        my_rank = df_sp100[df_sp100.record==str(record_time)].index[0] + 1
        
        return render_template('sp100.html', fam=fam, gen=gen, next_url=next_url, my_rank=my_rank, sp100=sp100)

    if (command == "stop"):
        sp100_mode = False

        print(sp100_mode)
        print(request.referrer)

        return redirect(request.referrer) 

@app.route('/stats')
def stats():
    global next_url
    
    fam = request.args.get('fam', default="", type=str)
    gen = request.args.get('gen', default="", type=str)
    next_url = request.args.get('next_url', default="", type=str)

    sql_str = 'SELECT timestamp, record FROM sp100 ORDER BY record'
    sp100 = query_db(sql_str)
        
    sql_str = 'SELECT timestamp, record FROM sp100 ORDER BY timestamp DESC'
    sp100_2 = query_db(sql_str)

    return render_template('stats.html', fam=fam, gen=gen, next_url=request.referrer, sp100=sp100, sp100_2=sp100_2)
"""

if __name__ == "__main__":
    # app.run(debug=True)
    app.run(debug=True, host='0.0.0.0')
