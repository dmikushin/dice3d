/**
 * @fileoverview Roll the die (3D version).
 * Assumes pre3d (http://deanm.github.com/pre3d/).
 * @author mimami24im@gmail.com
 */

/**
 * @namespace
 */
var midice3d = {};

/*constant*/
midice3d.con = {
  ttype : {
    // transform type
    t : 'translate', // mobile
    r : 'rotate',    // rotation
    s : 'scale'      // Expansion
  },
  rad1 : Math.PI / 180 // How many radians is 1 °
};

/** @typedef {{x: number, y: number, z: number}}*/
midice3d.Vector;
/**
 * Copy QuadFace
 * @param {Pre3d.QuadFace} qf Source QuadFace
 * @return {Pre3d.QuadFace} Copy result
 * @private
 */
midice3d.copyQF_ = function(qf) {
  var chkArr = [ qf.i0, qf.i1, qf.i2, qf.i3 ];
  var copiedArr = [];
  // i0, i1, i2, i3 is number for Shape.quads [],
  // Midice3d.Vector for Renderer.buffered_quads _ []. Qf
  for (var i = 0; i < chkArr.length; i++) {
    if (typeof chkArr[i] == 'object' && chkArr[i] !== null) {
      copiedArr[i] = midice3d.copyVec_(chkArr[i]);
    } else {
      copiedArr[i] = chkArr[i];
    }
  }
  var cqf = new Pre3d.QuadFace(copiedArr[0], copiedArr[1], copiedArr[2],
                               copiedArr[3]);
  cqf.centroid = (qf.centroid === null) ? null : midice3d.copyVec_(qf.centroid);
  cqf.normal1 = (qf.normal1 === null) ? null : midice3d.copyVec_(qf.normal1);
  cqf.normal2 = (qf.normal2 === null) ? null : midice3d.copyVec_(qf.normal2);
  return cqf;
};
/**
 * Copy Path
 * @param {Pre3d.Path} pt copy source path
 * @return {Pre3d.Path} Copy result
 * @private
 */
midice3d.copyPath_ = function(pt) {
  var cpt = new Pre3d.Path();
  for (var i = 0; i < pt.points.length; i++) {
    cpt.points[i] = midice3d.copyVec_(pt.points[i]);
  }
  for (var i = 0; i < pt.curves.length; i++) {
    cpt.curves[i] =
        new Pre3d.Curve(pt.curves[i].ep, pt.curves[i].c0, pt.curves[i].c1);
  }
  cpt.starting_point = pt.starting_point;
  return cpt;
};
/**
 * Copy Vector
 * @param {midice3d.Vector} vec Source Vector
 * @return {midice3d.Vector} Copy result
 * @private
 */
midice3d.copyVec_ = function(vec) { return {x : vec.x, y : vec.y, z : vec.z}; };
/**
 * Transform the value of Path. rotate rotates in the order of rotateZ → rotateY
 * → rotateX
 * @param {Pre3d.Path} path 対象path
 * @param {string} type The type of transform. Set the value of
 *     midice3d.con.ttype
 * @param {midice3d.Vector} tp transform size in each direction.
 * Rotation angle is specified in degrees
 * @private
 */
midice3d.transPath_ = function(path, type, tp) {
  var t = new Pre3d.Transform();
  switch (type) {
  case midice3d.con.ttype.t: {
    t.translate(tp.x, tp.y, tp.z);
    break;
  }
  case midice3d.con.ttype.r: {
    t.rotateZ(midice3d.con.rad1 * tp.z);
    t.rotateY(midice3d.con.rad1 * tp.y);
    t.rotateX(midice3d.con.rad1 * tp.x);
    break;
  }
  case midice3d.con.ttype.s: {
    t.scale(tp.x, tp.y, tp.z);
    break;
  }
  default: {
    return;
  }
  }
  if (path.starting_point === null) {
    var plen = path.points.length;
    path.points[plen] = {x : 0, y : 0, z : 0};
    path.starting_point = plen;
  }
  for (var i = 0; i < path.points.length; i++) {
    path.points[i] = t.transformPoint(path.points[i]);
  }
};

/**
 * Dice class.
 * Unless otherwise specified, the unit of length is a value where canvas.height
 * / 2 is 1.
 * @param {HTMLElement} canvas Display target HTMLCanvasElement
 * @param {number} width One side length
 * @constructor
 */
midice3d.Dice = function(canvas, width) {
  /**
   * One side length
   * @type {number}
   * @private
   */
  this.width_ = width;
  /**
   * Half the length of one side
   * @type {number}
   * @private
   */
  this.hfwdth_ = this.width_ / 2;
  /**
   * Display target HTMLCanvasElement
   * @type {HTMLElement}
   * @private
   */
  this.canvas_ = canvas;
  /**
   * Half of canvas.height (px)
   * @type {number}
   * @private
   */
  this.hfcvh_ = canvas.height / 2;
  /**
   * Eye colors other than 1. Set the character string including CSS color.
   * @type {string}
   * @private
   */
  this.eyecol_ = '#000000';
  /**
   * 1 eye color. Set the character string including CSS color.
   * @type {string}
   * @private
   */
  this.eye1col_ = '#ff0000';
  /**
   * An array of rotation amounts from the initial value of the object seen from
   * the camera. Create as many elements as there are buffers By default, the
   * camera looks down from above in the z-axis direction. Specify the angle of
   * rotation about each axis in degrees. rotateZ → rotateY → rotateX in order
   * @type {Array.<midice3d.Vector>}
   * @private
   */
  this.camAngles_ = [
    {x : 47, y : 40, z : 0}, {x : -40, y : 55, z : 0}, {x : 47, y : 35, z : -20}
  ];
  /**
   * The amount of movement of the object seen from the camera from the initial
   * value. By default, the camera looks down at {x: 0, y: 0, z: -1} directly
   * under the z axis. Specify the length to move in each axis direction.
   * @type {midice3d.Vector}
   * @private
   */
  this.camTrans_ = {x : 0, y : 0, z : -2};
  /**
   * Camera focal_length
   * @type {number}
   */
  this.camFocLen = 1;
  /**
   * Brightness contrast. Set a value between 0 and 1 inclusive.
   * The smaller the value, the higher the contrast.
   * @type {number}
   * @private
   */
  this.lightCont_ = 0.9;
  /**
   * The amount of rotation from the initial value of the light source. The
   * initial value is the unit vector in the z-axis direction. Specify the angle
   * of rotation about each axis in degrees. rotateX → rotateY → rotateZ in
   * order
   * @type {midice3d.Vector}
   * @private
   */
  this.lightAngle_ = {x : -30, y : 0, z : 20};
  /**
   * Light source direction vector. Generated by this.init ().
   * @type {midice3d.Vector}
   * @private
   */
  this.lightVec_ = {x : 0, y : 0, z : 1};
  /**
   * An array of dice display buffers. The size must be less than
   * this.camAngles_.length.
   * @type {Array.<midice3d.DBuffer_>}
   * @private
   */
  this.dbuf_ = [];
  /**
   * Renderer object. Generated by this.init ().
   * @type {Pre3d.Renderer}
   * @private
   */
  this.renderer_ = null;
};
/**
 * Initial processing. Create the dice display buffer in advance.
 */
midice3d.Dice.prototype.init = function() {
  // Renderer object generation
  var renderer = new Pre3d.Renderer(this.canvas_);
  renderer.camera.focal_length = this.camFocLen;
  renderer.fill_rgba = new Pre3d.RGBA(1, 1, 1, 1); // white

  // Create light source direction vector
  var g_z_axis_vector = {x : 0, y : 0, z : 1}; // z-axis direction unit vector
  var gzt = new Pre3d.Transform();
  gzt.reset();
  gzt.rotateX(midice3d.con.rad1 * this.lightAngle_.x);
  gzt.rotateY(midice3d.con.rad1 * this.lightAngle_.y);
  gzt.rotateZ(midice3d.con.rad1 * this.lightAngle_.z);
  this.lightVec_ = gzt.transformPoint(g_z_axis_vector);

  // Create a basic shape (die before rotation). Below, array index is (eye
  // value-1) Rotation angle to move from the z u003d 0 plane to the base plane
  var brot = [
    {x : -90, y : 0, z : 0}, {x : 0, y : -90, z : 0}, {x : 0, y : 0, z : 0},
    {x : 180, y : 0, z : 0}, {x : 0, y : 90, z : 0}, {x : 90, y : 0, z : 0}
  ];
  // Center of gravity of each side of basic form
  var hfwdth = this.hfwdth_;
  var cent = [
    {x : 0, y : hfwdth, z : 0}, {x : -hfwdth, y : 0, z : 0},
    {x : 0, y : 0, z : hfwdth}, {x : 0, y : 0, z : -hfwdth},
    {x : hfwdth, y : 0, z : 0}, {x : 0, y : -hfwdth, z : 0}
  ];
  // Path of the basic form of eyes 1 to 6
  var eyepathBase = [];
  var dice = new midice3d.DFace_(this.width_);
  for (var i = 0; i < 6; i++) {
    eyepathBase[i] = dice.makep(i + 1);
    for (var j = 0; j < eyepathBase[i].length; j++) {
      midice3d.transPath_(eyepathBase[i][j], midice3d.con.ttype.r, brot[i]);
    }
  }
  // Rotation angle required to raise the eye value
  var eyerot = [
    {x : 0, y : 0, z : 0}, {x : 0, y : 0, z : -90}, {x : -90, y : 0, z : 0},
    {x : 90, y : 0, z : 0}, {x : 0, y : 0, z : 90}, {x : 180, y : 0, z : 0}
  ];

  // Dice display buffer creation
  var cube, eyepath, pathidx, eyevals;
  for (var i = 0; i < this.camAngles_.length; i++) {
    renderer.emptyBuffer();
    var dicebuf = new midice3d.DBuffer_();
    // Camera settings according to the rotation angle of the dice
    this.set_camera_(renderer, this.camTrans_, this.camAngles_[i]);

    // Cube buffer creation
    cube = Pre3d.ShapeUtils.makeCube(this.width_ / 2);
    renderer.bufferShape(cube);
    this.chgIntensity_(renderer);
    dicebuf.drbuf = this.makeBuf_(renderer);

    // Create buffer for path with upper eyes 1 to 6
    for (var j = 0; j < 6; j++) {
      // Find the eye value to display
      eyevals = this.getDispEye_(cent, this.camAngles_[i], eyerot[j]);

      dicebuf.eyePathArr[j] = [];
      dicebuf.eye1idxArr[j] = [];
      pathidx = 0; // dicebuf.eyePathArr[j][] のindex
      // Creating a buffer for the displayed eye
      for (var j2 = 0; j2 < eyevals.length; j2++) {
        var eyeval = eyevals[j2];
        var pathlen = eyepathBase[eyeval].length;
        if (pathlen == 1) {
          // 1st eye index storage
          dicebuf.eye1idxArr[j].push(pathidx);
        }
        for (var j3 = 0; j3 < pathlen; j3++) {
          eyepath = midice3d.copyPath_(eyepathBase[eyeval][j3]);
          midice3d.transPath_(eyepath, midice3d.con.ttype.r, eyerot[j]);
          dicebuf.eyePathArr[j][pathidx] = this.makePBuf_(eyepath, renderer);
          pathidx++;
        }
      }
    }
    this.dbuf_[i] = dicebuf;
  }

  // Set Renderer object to property
  renderer.emptyBuffer();
  this.renderer_ = renderer;
};
/**
 * Calculate the eye value to display
 * @param {Array. <midice3d.Vector>} centroids centroid of face
 * @param {midice3d.Vector} camrot Camera rotation angle
 * @param {midice3d.Vector} eyerot Rotation angle to raise the eye value
 * @return {Array. <number>} eye value to display-1
 * @private
 */
midice3d.Dice.prototype.getDispEye_ = function(centroids, camrot, eyerot) {
  var eyeval = [];
  // Create transformation matrix
  var trans = new Pre3d.Transform;
  trans.rotateZ(midice3d.con.rad1 * eyerot.z);
  trans.rotateY(midice3d.con.rad1 * eyerot.y);
  trans.rotateX(midice3d.con.rad1 * eyerot.x);
  trans.rotateZ(midice3d.con.rad1 * camrot.z);
  trans.rotateY(midice3d.con.rad1 * camrot.y);
  trans.rotateX(midice3d.con.rad1 * camrot.x);

  // Apply transformation matrix to centroids and extract centroids with z> 0
  var tp;
  for (var i = 0; i < centroids.length; i++) {
    tp = trans.transformPoint(centroids[i]);
    // Consider floating point error
    if (tp.z > 0.0000001) {
      eyeval.push(i);
    }
  }
  return eyeval;
};

/**
 * Create draw buffer from renderer.buffered_quads_
 * @param {Pre3d.Renderer} renderer
 * @return {Array.<{qf: Pre3d.QuadFace, fillcolor: string}>}
 * @private
 */
midice3d.Dice.prototype.makeBuf_ = function(renderer) {
  var drbuf = [];
  var all_quads = renderer.buffered_quads_;
  all_quads.sort(midice3d.zSorter_);
  var obj, qf, frgba, iy, rgba, pushedp, iarr, nextj;
  for (var i = 0; i < all_quads.length; i++) {
    obj = all_quads[i];

    // qf made
    qf = obj.qf;
    renderer.projectQuadFaceToCanvasIP(qf);
    for (var j = 0; j < 4; j++) {
      if (j == 3) {
        nextj = 0;
      } else {
        nextj = j + 1;
      }
      pushedp = midice3d.pushPoints2d_(qf['i' + j], qf['i' + nextj]);
      qf['i' + j] = pushedp[0];
      qf['i' + nextj] = pushedp[1];
    }
    // Convert all coordinates to integers
    for (var j = 0; j < 4; j++) {
      qf['i' + j].x = ~~qf['i' + j].x;
      qf['i' + j].y = ~~qf['i' + j].y;
    }
    drbuf[i] = {};
    drbuf[i].qf = midice3d.copyQF_(qf);

    // create fillcolor
    frgba = obj.fill_rgba;
    iy = obj.intensity;
    rgba = [
      ~~(frgba.r * iy * 255), ~~(frgba.g * iy * 255), ~~(frgba.b * iy * 255),
      frgba.a
    ];
    drbuf[i].fillcolor = 'rgba(' + rgba.join(',') + ')';
  }
  return drbuf;
};
/**
 * Pre3d.Renderer.buffered_quads_ u200bu200bZ-coordinate sort function.
 * zSorter () from pre3d.js.
 * @param {Object} x comparison
 * @param {Object} y comparison
 * @return {number} Positive if x> y, 0 if x u003du003d y, negative otherwise
 * @private
 */
midice3d.zSorter_ = function(x,
                             y) { return x.qf.centroid.z - y.qf.centroid.z; };
/**
* Separate the distance between two points.

* @param {{x: number, y: number}} a coordinate

* @return {Array. <{x: number, y: number}>} Converted coordinates. [0] is a and
[1] is b

*/
midice3d.pushPoints2d_ = function(a, b) {
  var vec = Pre3d.Math.unitVector2d(Pre3d.Math.subPoints2d(b, a));
  if (isNaN(vec.x) || isNaN(vec.y)) {
    vec = {x : 0, y : 0};
  }
  var bpushed = Pre3d.Math.addPoints2d(b, vec);
  var apushed = Pre3d.Math.subPoints2d(a, vec);
  return [ apushed, bpushed ];
};

/*

* @param {Pre3d.Path} path Creation target
* @param {Pre3d.Renderer} renderer makes the environment
* @return {Pre3d.Path} Path that has the coordinates transformed and
projectPointToCanvas

*/
midice3d.Dice.prototype.makePBuf_ = function(path, renderer) {
  var trans = new Pre3d.Transform();
  trans.multTransform(renderer.camera.transform);
  trans.multTransform(renderer.transform);
  var tps = [];
  for (var i = 0; i < path.points.length; i++) {
    tps[i] = trans.transformPoint(path.points[i]);
  }

  var screen_points = renderer.projectPointsToCanvas(tps);
  var pathbuf = midice3d.copyPath_(path);
  pathbuf.points = screen_points;
  if (path.starting_point === null) {
    var start_point = renderer.projectPointToCanvas(
        trans.transformPoint({x : 0, y : 0, z : 0}));
    var spidx = pathbuf.points.length;
    pathbuf.points[spidx] = start_point;
    pathbuf.starting_point = spidx;
  } else {
    pathbuf.starting_point = path.starting_point;
  }

  // Convert coordinates to integers
  for (var i = 0; i < pathbuf.points.length; i++) {
    pathbuf.points[i].x = ~~pathbuf.points[i].x;
    pathbuf.points[i].y = ~~pathbuf.points[i].y;
    pathbuf.points[i].z = ~~pathbuf.points[i].z;
  }
  return pathbuf;
};

/**
 * Set camera. Move in order of rotateZ → rotateY → rotateX → translate
 * @param {Pre3d.Renderer} renderer Set target Renderer
 * @param {midice3d.Vector} trans
 * Length to move the object in each axis as seen from the camera
 * @param {midice3d.Vector} rot
 * A value that rotates the object around each axis as seen from the camera.
 * Unit is radian
 * @private
 */
midice3d.Dice.prototype.set_camera_ = function(renderer, trans, rot) {
  var ct = renderer.camera.transform;
  ct.reset();
  ct.rotateZ(midice3d.con.rad1 * rot.z);
  ct.rotateY(midice3d.con.rad1 * rot.y);
  ct.rotateX(midice3d.con.rad1 * rot.x);
  ct.translate(trans.x, trans.y, trans.z);
};
/**
 * Change brightness
 * @param {Pre3d.Renderer} renderer Renderer to change
 * @private
 */
midice3d.Dice.prototype.chgIntensity_ = function(renderer) {
  var intensity; // Surface brightness
  var n1;        // Face normal vector. Used for brightness calculation
  for (var i = 0; i < renderer.buffered_quads_.length; i++) {
    n1 = renderer.buffered_quads_[i].qf.normal1;
    intensity = Pre3d.Math.dotProduct3d(this.lightVec_, n1);
    intensity += (1 - intensity) * this.lightCont_;
    // Change brightness
    renderer.buffered_quads_[i].intensity = intensity;
  }
};
/**
 * Show dice
 * @param {number} bufidx index of the buffer to display
 * @param {number} upeyeval The eye value above
 * @param {number} x Center X coordinate (upper left is (0,0). Canvas right
 *     direction is positive)
 * @param {number} y center Y coordinate (upper left is (0,0). Canvas downward
 *     direction is positive)
 * @param {boolean} clrflg true to clear the Canvas before drawing. Otherwise
 *     false
 */
midice3d.Dice.prototype.draw = function(bufidx, upeyeval, tx, ty, clrflg) {
  if (clrflg === undefined) {
    clrflg = true;
  }

  var renderer = this.renderer_;
  var dicebuf = this.dbuf_[bufidx];

  if (clrflg) {
    renderer.clearBackground();
  }
  var ctx = renderer.ctx;
  // Draw a cube
  ctx.save();
  ctx.translate(tx, ty);
  var qf;
  for (var i = 0; i < dicebuf.drbuf.length; i++) {
    qf = dicebuf.drbuf[i].qf;
    ctx.beginPath();
    ctx.moveTo(qf.i0.x, qf.i0.y);
    ctx.lineTo(qf.i1.x, qf.i1.y);
    ctx.lineTo(qf.i2.x, qf.i2.y);
    ctx.lineTo(qf.i3.x, qf.i3.y);
    ctx.fillStyle = dicebuf.drbuf[i].fillcolor;
    ctx.fill();
  }

  // Draw eyes
  var eyePath, key, startpt;
  // Combine the index of the first eye into one string
  var eye1idxstr = '#' + dicebuf.eye1idxArr[upeyeval - 1].join('#') + '#';
  for (var i = 0, il = dicebuf.eyePathArr[upeyeval - 1].length; i < il; i++) {
    eyePath = dicebuf.eyePathArr[upeyeval - 1][i];
    // Determine if path is 1
    key = '#' + i + '#';
    if (eye1idxstr.indexOf(key) >= 0) {
      ctx.fillStyle = this.eye1col_;
    } else {
      ctx.fillStyle = this.eyecol_;
    }
    // draw the path
    ctx.beginPath();
    startpt = eyePath.starting_point;
    ctx.moveTo(eyePath.points[startpt].x, eyePath.points[startpt].y);
    var curves = eyePath.curves;
    for (var j = 0; j < curves.length; j++) {
      var curve = curves[j];
      if (curve.isQuadratic() === true) {
        var c0 = eyePath.points[curve.c0];
        var ep = eyePath.points[curve.ep];
        ctx.quadraticCurveTo(c0.x, c0.y, ep.x, ep.y);
      } else {
        var c0 = eyePath.points[curve.c0];
        var c1 = eyePath.points[curve.c1];
        var ep = eyePath.points[curve.ep];
        ctx.bezierCurveTo(c0.x, c0.y, c1.x, c1.y, ep.x, ep.y);
      }
    }
    ctx.fill();
  }
  ctx.restore();
};

/**
 * Dice display buffer class.
 * @private
 * @constructor
 */
midice3d.DBuffer_ = function() {
  /**
   * The result of projectQuadFaceToCanvasIP () of qf of buffered_quads_
   * u200bu200b[], Value to set in ctx.fillStyle.
   * @type {Array.<{qf: Pre3d.QuadFace, fillcolor: string}>}
   */
  this.drbuf = [];
  /**
   * Path array of eyes to display. index is the value of the upper eye -1.
   * @type {Array.<Array.<Pre3d.Path>>}
   */
  this.eyePathArr = [];
  /**
   * index of the first eye in eyePathArr
   * @type {Array.<Array.<number>>}
   */
  this.eye1idxArr = [];
};

/**
 * Dice surface class.
 * Unless otherwise specified, the unit of length is canvas.height / 2 with a
 * value of 1.
 * @param {number} width Side length
 * @private
 * @constructor
 */
midice3d.DFace_ = function(width) {
  /**
   * One side length
   * @type {number}
   * @private
   */
  this.width_ = width;
  /**
   * Half the length of one side
   * @type {number}
   * @private
   */
  this.hfwdth_ = this.width_ / 2;
  /**
   * What is the radius of the circle other than 1 times the side length?
   * @type {number}
   * @private
   */
  this.eyeRadRate_ = 0.1;
  /**
   * How many times the radius of the 1st eye circle is larger than the radius
   * of the 1st eye circle
   * @type {number}
   * @private
   */
  this.eye1RadRate_ = 1.5;
  /**
   * Percentage of the part where the center of the 4th circle is shifted from
   * the 1/4 of one side to the outside. Set a value between 0 and 1 inclusive.
   * @type {number}
   * @private
   */
  this.cdist_ = 0.1;
  /**
   * Distance from the center of the 1st eye to the center of the 4th eye circle
   * on the x coordinate
   * @type {number}
   * @private
   */
  this.distCir_ = this.width_ / 4 * (1 + this.cdist_);
  /**
   * Radius of eye circles other than 1
   * @type {number}
   * @private
   */
  this.eyeRad_ = this.width_ * this.eyeRadRate_;
  /**
   * Radius of the 1st circle
   * @type {number}
   * @private
   */
  this.eye1Rad_ = this.eyeRad_ * this.eye1RadRate_;
  /**
   * Coordinate array of eye centers.
   * The array index is 1: upper right 2: middle right 3: lower right 4: center
   * 5: upper left 6: middle left 7: lower left. The array element of the second
   * layer is [x coordinate, y coordinate]. Let the center of the dice be (0,0)
   * @type {Array.<?Array.<number>>}
   * @private
   */
  this.cxyArr_ = [
    null, [ this.distCir_, this.distCir_ ], [ this.distCir_, 0 ],
    [ this.distCir_, -this.distCir_ ], [ 0, 0 ],
    [ -this.distCir_, this.distCir_ ], [ -this.distCir_, 0 ],
    [ -this.distCir_, -this.distCir_ ]
  ];
};
/**
 * Create eye path on z u003d 0 plane
 * @param {number} value The eye value. An integer from 1 to 6
 * @return {Array. <Pre3d.Path>} created path
 */
midice3d.DFace_.prototype.makep = function(value) {
  // Argument array of this.makeEyeP_ for each eye
  var eyeposArr = [
    null, null, [ 5, 3 ], [ 1, 4, 7 ], [ 1, 3, 5, 7 ], [ 1, 3, 4, 5, 7 ],
    [ 1, 2, 3, 5, 6, 7 ]
  ];
  var paths = [];
  switch (value) {
  case 1: {
    var cp = Pre3d.PathUtils.makeCircle();
    midice3d.transPath_(cp, midice3d.con.ttype.t,
                        {x : -this.hfwdth_, y : 0, z : this.hfwdth_});
    var eye1R = this.eye1Rad_ / 0.5; // radius of makeCircle () is 0.5
    midice3d.transPath_(cp, midice3d.con.ttype.s,
                        {x : eye1R, y : eye1R, z : 1});
    paths[0] = cp;
    break;
  }
  case 2:
  case 3:
  case 4:
  case 5:
  case 6: {
    for (var i = 0; i < eyeposArr[value].length; i++) {
      paths[i] = this.makeEyeP_(eyeposArr[value][i]);
    }
    break;
  }
  default: {
    break;
  }
  }
  return paths;
};
/**
 * Create path for die rolls other than 1.
 * @param {number} pos 1: top right 2: middle right 3: bottom right 4: center 5:
 *     top left 6: center left 7: bottom left
 * @return {Pre3d.Path} created path
 * @private
 */
midice3d.DFace_.prototype.makeEyeP_ = function(pos) {
  var cp = Pre3d.PathUtils.makeCircle();
  pos = ~~pos;
  if (pos >= 1 && pos <= 7) {
    // Get coordinates of eye center
    var cx = this.cxyArr_[pos][0];
    var cy = this.cxyArr_[pos][1];
    // Eye path creation
    midice3d.transPath_(cp, midice3d.con.ttype.t,
                        {x : -this.hfwdth_, y : 0, z : this.hfwdth_});
    var eyeR = this.eyeRad_ / 0.5;
    midice3d.transPath_(cp, midice3d.con.ttype.s, {x : eyeR, y : eyeR, z : 1});
    midice3d.transPath_(cp, midice3d.con.ttype.t, {x : cx, y : cy, z : 0});
  }
  return cp;
};
