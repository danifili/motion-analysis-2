import React, {Component} from "react"
import Draggable, {DraggableCore} from 'react-draggable';

const ConnectingLine = ({ from, to }) => (
  <line
    x1={from.x}
    y1={from.y}
    x2={to.x}
    y2={to.y}
    stroke="rgb(200, 200, 200)"
    strokeDasharray="2,2"
    strokeWidth={1.5}
  />
);

class BezierCurve extends Component {
  constructor(props) {
    super(props);
    this.state = {
      startPoint: {x: 25, y: 25},
      controlPoint: {x: 50, y: 60},
      endPoint: {x: 75, y: 75}
    }
  }

  path() {
    return (
      <path
        d={`
          M ${[this.state.startPoint.x, this.state.startPoint.y]}
          Q ${[this.state.controlPoint.x, this.state.controlPoint.y]} ${[this.state.endPoint.x, this.state.endPoint.y]}
        `}
        fill="none"
        stroke="hotpink"
        strokeWidth={1.5}
      />
    );
  }

  renderNumber(num) {
    return parseFloat(Math.round( num * 10) / 10).toFixed(1)
  }

  renderPointCoords(id, point) {
    return id + ": " + this.renderNumber(point.x) + " " + this.renderNumber(point.y)
  }

  render() {
    const bounds = {
      left: 2,
      top: 2,
      right: 98,
      bottom: 98
    }
    return (
      <div style={{display:"inline"}}>
        <svg
        viewBox="0 0 100 100"
        style={{ 
          backgroundImage: 'url(' + this.props.img + ')',
          backgroundSize: '100% 100%' ,
        }}
        >
          {this.path()}
          <ConnectingLine from={this.state.startPoint} to={this.state.controlPoint}/>
          <ConnectingLine from={this.state.endPoint} to={this.state.controlPoint}/>
          <Draggable 
          defaultPosition={{x: this.state.startPoint.x, y: this.state.startPoint.y}}
          onDrag={(e, data) => this.setState({startPoint: {x: data.x, y: data.y}})}
          bounds={bounds}
          scale={5}
          >
            <ellipse
            cx={0} cy={0}
            rx={1.5} ry={1.5}
            />
          </Draggable>
          <Draggable 
          defaultPosition={{x: this.state.controlPoint.x, y: this.state.controlPoint.y}}
          onDrag={(e, data) => this.setState({controlPoint: {x: data.x, y: data.y}})}
          bounds={bounds}
          scale={5}
          >
            <ellipse
            cx={0} cy={0}
            rx={1.5} ry={1.5}
            fill={"blue"}
            />
          </Draggable>
          <Draggable 
          defaultPosition={{x: this.state.endPoint.x, y: this.state.endPoint.y}}
          onDrag={(e, data) => this.setState({endPoint: {x: data.x, y: data.y}})}
          bounds={bounds}
          scale={5}
          >
            <ellipse
            cx={0} cy={0}
            rx={1.5} ry={1.5}
            />
          </Draggable>  
        }
        </svg>
        <div>
          <ol>
            <li>{this.renderPointCoords("start point", this.state.startPoint)}</li>
            <li>{this.renderPointCoords("end point", this.state.endPoint)}</li>
            <li>{this.renderPointCoords("control point", this.state.controlPoint)}</li>
          </ol>
        </div>
      </div>
    )
  }
}

export default BezierCurve;