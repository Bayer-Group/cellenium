import Fly from './svg/fly.svg';

interface Icon {
  size: number;
}
function FlyIcon({ size }: Icon) {
  return <img src={Fly} className="h-auto" style={{ width: size }} alt="fly icon" />;
}

export default FlyIcon;
